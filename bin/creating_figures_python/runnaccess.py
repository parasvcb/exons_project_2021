import subprocess
import sys,os
if len(sys.argv)!=3:
    print ('please enter correct args 1inpdir 2outdir')
prog,inpdir,outdir=sys.argv
for i in os.listdir(inpdir):
    fname=os.path.join(inpdir,i)
    outname=os.path.join(outdir,i+'.rsa')
    if not os.path.isfile(outname):
        t1=os.path.basename(fname)
        print (t1)
        flag=t1.split('.')[0]
        res=subprocess.check_output(['naccess %s'%fname], shell=True, stderr=subprocess.STDOUT, stdout=FNULL)
        sys.exit()
        os.remove(flag+'.asa')
        os.remove(flag+'.log')
        subprocess.call(['mv %s %s'%(flag+'.rsa', outname)], shell=True)
