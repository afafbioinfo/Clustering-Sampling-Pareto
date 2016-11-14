import FileFunctions as FF
import os, subprocess
import conf as conf

def StructSampling(Pathconstraints, numberofsruct, Tmpr,Extension):
    dir = 'OutputSamples' + str(numberofsruct)
    FF.CreateFold(dir)
    for Path in Pathconstraints:
        for filename in FF.GetListFile(Path, Extension):
            Input = Path + "/" + filename
            output = dir + '/' + os.path.splitext(filename)[0]
            Command = 'RNAsubopt --noLP -p ' + str(numberofsruct) + ' -s -T ' + str(Tmpr)
            if Path ==conf.PathConstrainteFile:
                Command += ' -C'
            if Path == conf.PathConstrainteFileShape:
                ShapeFile = conf.PathConstrainteFileShape + "/" + os.path.splitext(filename)[0] + 'Shape.txt'
                Command += ' --shape ' + ShapeFile
            subprocess.call(Command, stdin=open(Input, 'r'), stdout=open(output, 'wb'), shell=True)

    return dir



