#example run :
#Rscript UpdateGraftPoint.r treeGraftoutput2.txt specieslist.csv|specieslist.tsv [saveFolder|savelist.out|savelist.txt]
#can have either names of species list or name of files + name of species list
#no column names

#   Rscript UpdateGraftPoint.r loxaf-corrected.out loxaf-spec-corrected.csv

#output text file at the end. 
#get an number for new graft point

library(stringr)
args = commandArgs(trailingOnly=TRUE)

args = c('twoToTest.csv','rainbowandfrogSpeciesOnly.csv') 

##################################################
############# verify/prepare inputs ##############
##################################################
if(length(args) == 3){
    
    if(grepl('/$',args[3])) {
        saveloc = args[3]
    } else {
        saveloc = str_c(args[3],'/')
    }
    tryCatch(
        {
            write.csv(data.frame(),str_c(saveloc,'updatedTreeGrafterPoint.out'))
        },
        error=function(e) {
            message('Save location invalid')
            quit()
        },
        warning=function(w) {
            message('Save location warning')
            print(w)
        }
    )     
} else if (length(args) != 2) {
    print('invalid number of inputs, try : Rscript UpdateGraftPoint.r treeGraftoutput.txt specieslist.csv [saveloc|savelist.out|savelist.txt] ')
    quit()
} else{
    saveloc = './updatedTreeGrafterPoint.out'
}

#make sure the inputs are correct
if(!file.exists(args[2])){
    print('Species input file does not exist. Quitting run. ')
    quit()
}
if(grepl('.csv$',args[2])) {
    inputDF = read.csv(args[2],header=F)
    } else if(grepl('.tsv$',args[2])) inputDF = read.csv(args[2],sep='\t',header=F)
namedCols = FALSE
if(dim(inputDF)[2] == 2)  {
    inputSpecies = inputDF[,2]
    namedCols = TRUE
    } else if(dim(inputDF)[2] == 1) {
    inputSpecies =  inputDF[,1] 
    } else {
    print('Species Input file has incorrect number of columns. Quitting run. ')
    quit()
}

if(file.exists(args[1])) {
    treeGrafterRes = readLines(args[1])   
    treeGrafterRes = treeGrafterRes[which(treeGrafterRes!="")] 
} else{
    print('output results file not found. Quitting run.')
    quit()
}

if(length(treeGrafterRes) > length(inputSpecies)){
    #print(length(treeGrafterRes))
    #print(length(inputSpecies))
    print('Length of species input does not match length of treeGrafter input. Quitting run. ')
    quit()
}
if(length(inputSpecies) == 0){
    print('No input species in file? Quitting run.')
    quit()
}


#getting list of AN numbers and pthr trees from the treeGrafter results
tgr = str_split(treeGrafterRes,' +|\t|,',simplify=T) # dataframe, name, pthr tree, pthrtree:SF# , go terms
tgr = gsub('"','',tgr)
print(dim(tgr))
print(head(tgr))
pthr = tgr[,2]
an = tgr[,3]
sampleName = tgr[,1]
unqpthr = unique(pthr) #don't want to reload pthr trees multiple times, so check for each unique panther tree
idNumber = 1:length(an) #used to keep records in order at end

if(sum(grepl('AN\\d+',an) == T) != length(an)){
    print('AN list in output file contains non AN number(s)')
    print(an[!grep('AN\\d+',an)])
    quit()
}
if(sum(grepl('PTHR\\d\\d\\d\\d\\d',pthr) == T) != length(pthr)){
    print('PTHR list in output file contains non PTHR number(s)')
    quit()
}

#TODO fix this, it doesn't work correctly
if(namedCols){
   mnc = match(inputDF[,1],sampleName)
   if(sum(is.na(mnc)) > 0){
       print('Some sample names in species list do not match sample names in tree grafter output list.')
       #quit()
   } else{
       inputSpecies = inputSpecies[mnc[!is.na(mnc)]]
   }
}

##################################################
########### Load in Needed Datasets ##############
##################################################

#child parent has the ncbi Taxon id numbers for the children and the closest ancestors in the panther taxon tree
childparent = read.csv('Data/NCBIChildParentOnly.csv')[,-1]

#short names is the 5 letter species designations
shortnames  = read.csv('Data/ncbitaxonshortnames.csv')[,-1] #remove column of column numbers

#speciesNamesToID gives the taxon id numbers for the species names in the panther taxon tree
speciesNamesToID = read.csv('Data/pthrNamesToTaxonIDs.csv')[,-1] 

##################################################
#Convert Input Species to Closest parent taxon ID#
##################################################

pnames = c() #pnames are the ncbitaxonids that are closest parent (or self) to input species, 'error' if 1 or id not recognized
for(i in 1:length(inputSpecies)){
    iname = inputSpecies[i]
    if(iname %in% shortnames[,1]) iname = shortnames[match(iname,shortnames[,1]),2]
    m = match(iname,childparent[,1])
    if(is.na(m)) {
        print(paste('ERROR, NCBITaxon ID not recognized for input species : ',iname, ' - Will be labeled as error in results with no graph output'))
        pname = 'error'    
    } else{
        pname = childparent[m,2] #closest ancestor in panther taxon tree to input species
        if(pname == 1) { #if parent of the input species is 1, output error message
            pname = 'error'
            print(paste('LUCA was the closet ancestor found for input species : ',iname, ' - in taxon tree, this species may be a virus, results omitted'))
        }
    }
    pnames = c(pnames,pname)
}

##################################################
########### Update The Graft Point ###############
##################################################

#print(sampleName)

newGraft = c() #list of the new nodes taken for the graft point
pathways = rep(c(),length(pthr)) #list the pathway between the original graft point and the new graft point for each input
orthos = rep(c(),length(pthr)) #the list of uniprot ids for orthologs for each input
orthosPaths  = rep(c(),length(pthr)) #the list of edges leading to orthologs for each input
newAN = rep('',length(pthr))
distFloat = rep('',length(pthr))
distBranches = rep('',length(pthr))

print(unqpthr)

source('treeGrafterUpgradeFunctions.r') #contains the read_panther function to read the tree with additional information, updated from the aphylo package to correct several bugs

for(i in 1:length(unique(unqpthr))){ 
    graftAN = an[pthr == unqpthr[i]] #this is AN number of treeGraft graft point, graftID for taxonID assigned later
    sampNames = sampleName[pthr == unqpthr[i]] #the assigned name from output file
    idNum = idNumber[pthr == unqpthr[i]] #used to keep track of the records to output in order input
    newGraftID = pnames[pthr == unqpthr[i]] 
    if(file.exists(str_c('PANTHER17.0_data/Tree_MSF/',unqpthr[i],'.tree'))){
        tree = read_panther(str_c('PANTHER17.0_data/Tree_MSF/',unqpthr[i],'.tree')) #function loaded from source file
    } else {
        print(paste('Panther Tree not found??? : ',unqpthr[i]))
        next
    }
    uniprot = str_split(tree$tree$tip.label,'UniProtKB=',simplify=T)[,2]

    leaves = str_split(tree$tree$tip.label,':',simplify=T)
    anlist = leaves[,1] #AN numbers of leaves
    leaves = str_split(leaves[,2],'\\|',simplify=T) 
    specleaves = leaves[,1] #gives 5 letter species name
    specInternal = tree$internal_nodes_annotations$ancestor #gives full species name
    specIDleaves = speciesNamesToID[match(specleaves,speciesNamesToID[,2]),3] #convert name to ID number
    specIDinternal = speciesNamesToID[match(specInternal,speciesNamesToID[,1]),3] 
    #note some internal are NA if they are duplications and thus don't have a species
    specIDinternal[is.na(specIDinternal)] = '-1'
    numLeaves = length(specleaves)

    countLeaves = 1:numLeaves
    countInternal = 1:length(specInternal)

    #creating secondary tree with length 1 edges, used to determine the distance between any two nodes, which will be initial method of determining closest new graft point if there are duplicate leaves
    tree2 = tree
    tree2$tree$edge.length = rep(1,length(tree$tree$edge.length))
    treeD = dist.nodes(tree$tree)
    tree2D = dist.nodes(tree2$tree)

    tree$tree$tip.label = uniprot #changed so labels on graphs just show uniprot id

    for(j in 1:length(graftAN)){
        print(paste(i,j))
        #print(sampNames[idNum[j]])
         
        #get graft taxon ID from an number, and get it's position in full node list to be used in the distance matrix
        graftNode = -1
        if(graftAN[j] %in% rownames(tree$internal_nodes_annotations)) { #internal node case
            graftID = speciesNamesToID[match(specInternal[match(graftAN[j],rownames(tree$internal_nodes_annotations))],speciesNamesToID[,1]),3] #get the id number for the specific graft point
            #print(paste('a',i,j))
            if(is.na(graftID)){ #for dealing with treegrafter Assigning to duplication nodes
                graftID = '-1'
                if(FALSE){ #this code is for the case where we want to update treegrafters initial node from a duplication, but this seems unnecessary
                    mgraft = match(graftAN[j],rownames(tree$internal_nodes_annotations))
                    count = 0
                    flag = FALSE
                    while(TRUE){
                        print(mgraft)
                        count = count + 1
                        mgraftPos = mgraft + numLeaves
                        mgraft = tree$tree$edge[match(mgraftPos,tree$tree$edge[,2]),1] - numLeaves
                        if(!is.na(tree$internal_nodes_annotations$type[mgraft] ) && tree$internal_nodes_annotations$type[mgraft] == 'S') {
                            graftID = speciesNamesToID[match(tree$internal_nodes_annotations$ancestor[mgraft],speciesNamesToID[,1]),3] 
                            break
                        }
                        if(count > 30 || mgraft == 1) {
                            
                            flag2 = TRUE
                            break
                        }
                    }
                    if(flag2) {
                        print('Infinite loop dealing with graft point assigned to duplication seemingly or duplication nodes to top')
                        newGraft[idNum[j]] = 'error'
                        next
                    }
                }
            }
            graftIDPos = numLeaves + match(graftAN[j],rownames(tree$internal_nodes_annotations))
        } else if(graftAN[j] %in% anlist) { #leaf node case
            #print(paste('b',i,j,specleaves[match(graftAN[j],anlist)]))
            #print(head(speciesNamesToID))
            #print(dim(speciesNamesToID))
            #print(graftAN[j] %in% anlist)
            #print(match(graftAN[j],anlist))
            #print(speciesNamesToID[111,])
            graftIDPos = match(graftAN[j],anlist)
            graftID = speciesNamesToID[match(specleaves[graftIDPos],speciesNamesToID[,2]),3]
            #print(paste(graftIDPos,graftID))
            } else { #not in leaf or internal node????
            print(paste('Error, no taxon id found for initial graft point for : ',sampNames[j]))
            newGraft[idNum[j]] = 'error'
            next
        }
        #print(graftID)
        #print(newGraftID)
        #Comparing the treeGrafter graft node with the closest graph node to input species
        if(graftID == newGraftID[j]){ #orthos still assigned so no 'next' command
            newGraft[idNum[j]] =  'identical'
            newAN[idNum[j]] = graftAN[j]
            pathways[[idNum[j]]] = NA
            currLoc = graftIDPos
            graftNode = graftIDPos
        } else if( newGraftID[j] == 'error') {
            #print('Propagating error, skipping')
            newGraft[idNum[j]] = 'error'
            pathways[[idNum[j]]] = NA
            orthos[[idNum[j]]] = uniprot
            orthosPaths[[idNum[j]]] = NA
            graftNode = length(spec) + 1 #TODO check this is right
            next
        } else{ #getting closest node that matches newGraftID[j] based on the distance based on branch numbers. Branch length used if duplicates.
            currGraftID = newGraftID[j]
            if(!(currGraftID %in% specIDleaves) && !(currGraftID %in% specIDinternal)) {
                print(paste('Input Species not in internal or external node species : ',sampNames[j]))
                count = 0
                flag = FALSE
                currGraftIDparentSpec = currGraftID
                while(TRUE){ 
                    count = count + 1
                    if(count == 30) break
                    #print(currGraftIDparentSpec)
                    testna = speciesNamesToID[match(currGraftIDparentSpec,speciesNamesToID[,3]),4] 
                    if(!is.na(testna)) {
                        currGraftIDparentSpec = testna
                    } else{
                        break
                    } 
                    #if(i == 5) print(speciesNamesToID)
                    #print(NA %in% tree$internal_nodes_annotations$ancestor)
                    #print(unqpthr[i])
                    if(!is.na(currGraftIDparentSpec) && currGraftIDparentSpec %in% tree$internal_nodes_annotations$ancestor) {
                        print(paste(currGraftIDparentSpec,'option1'))
                        break
                    }
                    if( currGraftIDparentSpec == specIDinternal[1]){ #if root is closest node, it would seemingly have no orthologs, only children
                        newGraft[idNum[j]] = 'root'
                        pathways[[idNum[j]]] = NA
                        orthos[[idNum[j]]] = uniprot
                        orthosPaths[[idNum[j]]] = NA
                        flag = TRUE
                        break
                    } 
                }
                if(flag) { 
                    print(paste('Root was the closet ancestor found for input species in taxon tree for sample : ',sampNames[j], ' - No output generated.'))
                    next
                }
                 currGraftID = currGraftIDparentSpec
            }
            #print(currGraftID)
            #if(i == 5) print(currGraftIDparentSpec) 
            # Determining where the matching nodes are located (using total node number (for internal nodes that total # leaves + internal number)) 
            #print(specIDleaves)
            ##print(currGraftID)
            #print(specIDleaves == currGraftID)
            matchGraft = countLeaves[specIDleaves == currGraftID]  
           
            if(length(matchGraft) == 0){
                matchGraft = countInternal[specIDinternal == currGraftID]
                if(length(matchGraft)==0) {
                    print("Error, no match to leaf or internal nodes found for new graft")
                    newGraft[idNum[j]] = 'error'
                    next
                }
                matchGraft = matchGraft + numLeaves
            }
            #print(matchGraft)
            #getting distance between graftIDPos (from tree grafter) and distance to all possibilities
            dbranchNums = tree2D[graftIDPos,matchGraft]
            #print(length(dbranchNums))
            #finding if more than one minimium distance and determining which has shortest overall branch length if so
            mind = min(dbranchNums)
            possibles = c()
            #if(i == 5) print(paste(k, dim(tree2D)[1],dim(tree2D)[2],graftIDPos,matchGraft[k],length(matchGraft),matchGraft))
            if(i == 5) print(paste(i,j,length(matchGraft)))
            for(k in 1:length(matchGraft)){ 
                if(tree2D[graftIDPos,matchGraft[k]] == mind) {
                    possibles = c(possibles,k)         
                }
            }
            distBranches[idNum[j]] = mind
            #print('possibiles')
            #print(possibles)
            if(length(possibles)==1){
                graftNode = matchGraft[possibles[1]]
                distFloat[idNum[j]] = treeD[graftIDPos,matchGraft[possibles[1]]]
            } else if(length(possibles) > 1){
                dbranchLens = treeD[graftIDPos,matchGraft[possibles]]
                mind = min(dbranchLens)
                distFloat[idNum[j]] = mind
                for(k in 1:length(possibles)){
                    if(treeD[graftIDPos,matchGraft[possibles[k]]] == mind) {
                        graftNode = matchGraft[possibles[k]]      
                    }
                }
            } else {
                    print("Error, 0 possibilities for minimum distance, ???")
                    newGraft[idNum[j]] = 'error'
                    next
            }


            #get pathway between the two graft points
            currLoc = graftNode #note graftNode is position in node list of updated graft point
            pathway = c()
            graftCurr = graftIDPos

            graftIDPosParents = c(graftCurr)
            graftParentsPath = c( )
            while(graftCurr != (numLeaves + 1)) { #get list of parents of graft point
                graftParentsPath = c(graftParentsPath,match(graftCurr,tree$tree$edge[,2]))
                graftIDPosParents = c(graftIDPosParents,graftCurr)
                graftCurr = tree$tree$edge[match(graftCurr,tree$tree$edge[,2]),1]
            }
            #print(paste(currLoc,graftNode))
            #print(numLeaves)
            while(currLoc != (numLeaves + 1)) { #get parents of input species until it matches a graft id parent
                mEdge = match(currLoc,tree$tree$edge[,2])
                currLoc = tree$tree$edge[mEdge,1]
                #if(currLoc %in% graftIDPosParents){
                #    break
                #}
                pathway = c(pathway,mEdge)
            }
            pathwaycombine = c(pathway,graftParentsPath) #remove the duplicates as those are overlapping lines
            pathwaycombine = as.integer(names(table(pathwaycombine))[table(pathwaycombine)==1])
            #get the orthologs of all spec parents
            newGraft[idNum[j]] = graftNode
            if(graftNode > numLeaves) {
                newAN[idNum[j]] = rownames(tree$internal_nodes_annotations)[graftNode - numLeaves]
            } else {
                newAN[idNum[j]] = anlist[graftNode]
            }
            pathways[[idNum[j]]] = pathwaycombine

        }
        
        ######### fix this later to add the children
        ### add all children nodes if graft point isn't same as input : 
        #if(newGraft[idNum[j]] != 'identical'){
         #   extract.clade(, newAN[idNum[j]])
        #}

        #want to find speciation nodes : 
        currLocParent = graftNode
        currLocParents = c()
        #currLocParentsPath = c( )
        while(currLocParent != (numLeaves + 1)) { 
            #currLocParentsPath = c(currLocParentsPath,match(currLocParent,tree$tree$edge[,2]))    
            currLocParent = tree$tree$edge[match(currLocParent,tree$tree$edge[,2]),1]
            currLocParents = c(currLocParents,currLocParent)
        }

        #want the parents that are speciation nodes so that we can get orthogonal proteins (not from a duplication)
        if(FALSE){

        validParents = currLocParents[tree$internal_nodes_annotations$type[currLocParents - numLeaves] == 'S']
        #make sure paths are saving all the correct edges
        if(length(validParents) == 0){ 
            print(paste('No Valid Parents that are speciation nodes? : ',sampNames[j]))
            next
        }

        #want to find all orthologs that first join up at a speciation node
        allPathsLeadToOrthologs = c()
        validLeaves = c()
        
        currLocParents = c(currLocParents,graftNode) #test
        print(graftNode)
        for(l in 1:length(anlist)){ #might be a more efficent way to do this, but duplication nodes create an issue, so will check all leaves
        #alternate method might be to get subtree for node + all duplications + highest parent, take all from parent and subtract other two
            leafLineLoc = l
            leafPath = c()
            while(leafLineLoc != (numLeaves + 1)) {  
                leafPath = c(leafPath,match(leafLineLoc,tree$tree$edge[,2]))    
                leafLineLoc = tree$tree$edge[match(leafLineLoc,tree$tree$edge[,2]),1]
                #test removing this and adding graftnode to currLocParents if(leafLineLoc == graftNode) break
                if(leafLineLoc %in% currLocParents){
                    if(leafLineLoc %in% validParents){
                        allPathsLeadToOrthologs = c(allPathsLeadToOrthologs,leafPath)
                        validLeaves = c(validLeaves,l)
                    }
                    else {
                        break
                    }
                }
            }
        }
        }

        if(TRUE){
            validParents = currLocParents[tree$internal_nodes_annotations$type[currLocParents - numLeaves] != 'S']
            #make sure paths are saving all the correct edges
            if(length(validParents) == length(currLocParents)){ 
                print(paste('No Valid Parents that are speciation nodes? : ',sampNames[j]))
                next
            }

            allPathsLeadToOrthologs = c()
            invalidLeaves = c()
            

            for(l in 1:length(anlist)){ #might be a more efficent way to do this, but duplication nodes create an issue, so will check all leaves
            #alternate method might be to get subtree for node + all duplications + highest parent, take all from parent and subtract other two
                leafLineLoc = l
                leafPath = c()
                while(leafLineLoc != (numLeaves + 1)) {  
                    leafPath = c(leafPath,match(leafLineLoc,tree$tree$edge[,2]))    
                    leafLineLoc = tree$tree$edge[match(leafLineLoc,tree$tree$edge[,2]),1]
                    #test removing this and adding graftnode to currLocParents if(leafLineLoc == graftNode) break
                    if(leafLineLoc %in% currLocParents){
                        if(leafLineLoc %in% validParents){
                            allPathsLeadToOrthologs = c(allPathsLeadToOrthologs,leafPath)
                            invalidLeaves = c(invalidLeaves,l)
                        }
                        else {
                            break
                        }
                    }
                }
            } 

            if(length(invalidLeaves) > 0) validLeaves = seq(1,length(anlist),1)[-invalidLeaves]
            else validLeaves = 1:length(anlist)
        }
        #Update initial lists, dealing with case for graphing in which the two graft points are the same
        validLeaves = unique(validLeaves)

        allPathsLeadToOrthologs = unique(allPathsLeadToOrthologs)
        if(newGraft[idNum[j]] ==  'identical') { #if new graft point is the same as old, we don't want first branch in path or would highlight itself as ortholog in graph
            matchIdent = match(match(graftIDPos,tree$tree$edge[,2]), allPathsLeadToOrthologs)
            if(!is.na(matchIdent)){
                allPathsLeadToOrthologs = allPathsLeadToOrthologs[-matchIdent]
            }
            mleaves = match(graftIDPos,validLeaves)
            if(!is.na(mleaves)){    
                validLeaves = validLeaves[-mleaves]
            }
        }

        #get the children of the graftIDPos
        #print('a')
        if(FALSE){
        toCheckList = c()
        checkCounter1 = 1
        currPos = graftIDPos
        childLeaves = c()
        while(TRUE){
            if(sum(tree$tree$edge[,1] == currPos) > 0)  {
                toCheckList = c(toCheckList,tree$tree$edge[tree$tree$edge[,1] == currPos,2])
            } else {
                childLeaves = c(childLeaves,currPos)
            }
            if(checkCounter1 == length(toCheckList) || length(toCheckList) == 0) break

            currPos = toCheckList[checkCounter1]
            checkCounter1 = checkCounter1 + 1
            #print(length(toCheckList))
            #print(length(unique(toCheckList)))
            #print(checkCounter1)

        }
        validLeaves = c(validLeaves,childLeaves)
        }
        validSpec = specleaves[validLeaves]
        validUni =  uniprot[validLeaves]
        validUni = validUni[!(validSpec == inputSpecies[idNum[j]])]
        #print('b')
        


        orthos[[idNum[j]]] = validUni
        orthosPaths[[idNum[j]]] = unique(allPathsLeadToOrthologs)

        ##################################################
        ################# Create Plot ####################
        ##################################################
        edgeColor = rep('#000000',length(tree$tree$edge))
        #edgeGreen = #unique(orthosPaths[[idNum[j]]])
        edgeOrange = unique(pathways[[idNum[j]]])
        #edgeColor[edgeGreen] = '#02FC0F'
        edgeColor[edgeOrange] = '#F98900'
        #print(head(edgeColor,30))
        nodeColor = rep('#000000',length(tree$tree$node.label))
        nodeColor[tree$internal_nodes_annotations$type != 'S'] = '#FC02FC'
        leafColor = rep('#000000',length(tree$tree$tip.label))
        leafColor[validLeaves] = '#02fc0f'
        leafColor[validLeaves & (specleaves[validLeaves] ==  inputSpecies[idNum[j]])] = '#FF0000'
        if(newGraft[idNum[j]] ==  'identical') leafColor[graftIDPos] = '#F98900' 
        pdf(str_c(saveloc,unqpthr[i],'-',graftAN[j],'-',inputSpecies[idNum[j]],'.pdf'),height = 1/12 * (length(tree$tree$tip.label)+length(tree$tree$node.label)) + 5)
            par(mar = c(c(8, 3, 3, 3) ))
            plot.phylo(tree$tree, main = str_c('Graft of ',sampNames[j],'\n into Panther Tree ',unqpthr[i]),edge.color = edgeColor, tip.color = leafColor, label.offset = 0.2,show.node.label=T,show.tip.label=T,align.tip.label=T) 
            nodelabels(pch = 21, bg = nodeColor, cex=1)
            legend(x=0,y=0,xpd = TRUE, lwd = 2, lty = c(NA,NA,1,1), bty='n', legend=c('Duplication','Speciation','Graft Point Movement','Ortholog Paths'), pch=c(21,21,NA,NA),pt.bg = c('#FC02FC','#000000','#F98900','#02FC0F'),col = c('#FC02FC','#000000','#F98900','#02FC0F'))
        dev.off()
        #print(paste(i,j,'done'))
    }
}

#print(data.frame(sampleName,newAN,newGraft))
uniprotList = rep('',length(pthr))
for(i in 1:length(orthos)){
    uniprotList[i] = str_c(orthos[[i]],collapse=';')
}
results = data.frame(sampleName,pthr,an,newAN,newGraft,distFloat,distBranches,uniprotList)
write.table(results,'./twoSpeciesTest.out',sep='\t')#str_c(saveloc,'loxaf-orthos-1500.out'),sep='\t')

