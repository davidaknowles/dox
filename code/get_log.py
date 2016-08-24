import sys
import glob
import os

basedir="/afs/cs.stanford.edu/u/davidknowles/scailscratch/dox/fastq/"

outdirs=glob.glob( basedir+"*_out" )

first=True

#fout=open("log.csv","w")

sys.stderr.write( "Num dirs: %d \n" % len(outdirs) )

for outdir in outdirs:
    fn=outdir + "/Log.final.out"
    if not os.path.isfile(fn):
        sys.stderr.write( "No "+fn+"\n")
        continue
    with open( fn, "r" ) as f:
        names=[]
        values=[]
        for l in f:
            if not "|" in l:
                continue
            l=l.split("|")
            n=l[-1].strip()
            n=n.strip("%")
            if len( n.split() )>1:
                continue
            values.append( float(n) )
            names.append( l[0].strip() )
        if first:
            print( "Sample\t" + "\t".join(names) )
            first=False
        print( outdir.split("/")[-1][0:-4] + "\t" + "\t".join( [str(g) for g in values] ) )
            
#fout.close()
