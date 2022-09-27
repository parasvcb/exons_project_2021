import constructing_data.mappings_refseqProtein as map
import sys,os
if len (sys.argv)!=6:
    print ('Please enter correct args 1. refseqDir, 2. prevSSDir 3. storeAssignDir 4. unassignedPendingDir 5. mapfilename')
prog,refDir,prevss,curss,pendingDir,fname=sys.argv


#print (refDir,prevss,curss,pendingDir,fname)
map.mapping_alreadyDone_vs_Pending(unassignedDir=refDir,alreadydoneDir=prevss,copyDonetoUnassigned=curss,pendingToRun=pendingDir,createTheMapFile=fname)