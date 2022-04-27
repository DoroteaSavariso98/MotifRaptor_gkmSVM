import sys
import argparse
import os
import subprocess
import time, datetime
import glob
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description='This script takes as input a vcf file of SNPs (or an input_snp.tsv in format chrX_POS_REF_ALT if no input is given) and returns the deltaSVM scores for each SNP calculated with the models of the TFs given in input (or all of them if no input is given)')
    parser.add_argument("-b", "--background", dest="background", required=False, help="vcf file of whole genome SNPs")
    parser.add_argument("-wd", "--wd", dest="wd", required=True, help="gkm_score directory")
    parser.add_argument("-m", "--models", dest="models", required=False, help="names of the TFs comma-separated")
    parser.add_argument("-o", "--outdir", dest="outdir", required=True, help="output directory")
    args = parser.parse_args()

    outdir = args.outdir
    gkm_folder = args.wd

    if args.background != None:
        with open (args.background, "r") as file:
            lines = file.readlines()
            file3 = open(gkm_folder+'/input_snp.tsv', "w")
            for line in lines:
                if line.startswith('#'):
                    continue 
                chrom, pos, snpID, ref, alt = line.rstrip().split("\t")[0:5]
                file3.write("{0}_chr{1}_{2}_{3}_{4}\n".format(snpID, chrom, pos, ref, alt))
            file3.close()


    tfs = []
    if args.models=="all":
        for i in os.listdir(gkm_folder+"/gkmsvm_models"):
            if i.endswith(".model.txt"):
                tfs.append(i.rstrip().split("_")[0])
    else:
        for i in (args.models).rstrip().split(","):
            tfs.append(i.upper())


    print("Command: generating score files for SNPs...")

    subprocess.call('set -e', shell=True)
    subprocess.call('rm -rf '+gkm_folder+'/data ' +outdir+ '/motifscanfiles', shell=True)
    subprocess.call('mkdir '+gkm_folder+'/data ' +outdir+ '/motifscanfiles', shell=True)
    print("Generating allelic sequences...")
    subprocess.call('python3 '+gkm_folder+'/scripts/generate_allelic_seqs.py -f '+gkm_folder+'/genome/hg19.fa -s '+gkm_folder+'/input_snp.tsv -o '+gkm_folder+'/data/selex_allelic_oligos', shell=True)
    df = pd.DataFrame()
    for tf in tfs:
        subprocess.call('mkdir '+gkm_folder+'/out '+gkm_folder+'/log', shell=True)
        print("Generating binding scores for {0}...".format(tf))
        model = os.path.basename(''.join(glob.glob(gkm_folder+'/gkmsvm_models/' + tf + '*.model.txt')))
        subprocess.call(gkm_folder+'/scripts/gkmpredict '+gkm_folder+'/data/selex_allelic_oligos.ref.fa '+gkm_folder+'/gkmsvm_models/' + model + ' '+gkm_folder+'/out/' + tf + '.ref.gkm.tsv &>'+gkm_folder+'/log/gkmpredict.ref.log', shell=True)
        subprocess.call(gkm_folder+'/scripts/gkmpredict '+gkm_folder+'/data/selex_allelic_oligos.alt.fa '+gkm_folder+'/gkmsvm_models/' + model + ' '+gkm_folder+'/out/' + tf + '.alt.gkm.tsv &>'+gkm_folder+'/log/gkmpredict.ref.log', shell=True)
        subprocess.call('paste '+gkm_folder+'/out/' + tf + '.ref.gkm.tsv '+gkm_folder+'/out/' + tf + '.alt.gkm.tsv > '+gkm_folder+'/out/' + tf + '.gkm.tsv', shell=True)  
        file = pd.read_csv(gkm_folder+'/out/' + tf + '.gkm.tsv', sep='\t', header=None)
        df['snpID'] = file[0]
        df['binding_ref'] = file[1]
        df['binding_alt'] = file[3]
        print("Generating deltasvm scores for {0}...".format(tf))
        weightfilename = os.path.basename(''.join(glob.glob(''+gkm_folder+'/weights/*' + tf + '.weights.txt')))
        subprocess.call('python3 '+gkm_folder+'/scripts/deltasvm.py -r '+gkm_folder+'/data/selex_allelic_oligos.ref.fa -a '+gkm_folder+'/data/selex_allelic_oligos.alt.fa -w '+gkm_folder+'/weights/' + weightfilename + ' -o '+gkm_folder+'/out/' + tf + '.delta.tsv &>'+gkm_folder+'/log/deltasvm.log', shell=True)
        file2 = pd.read_csv(gkm_folder+'/out/' + tf + '.delta.tsv', sep='\t', header=None)
        df['deltaSVM_score'] = file2[1]
        subprocess.call('rm -rf '+gkm_folder+'/out '+gkm_folder+'/log', shell=True)
    df.to_csv(outdir+ '/motifscanfiles/' + tf + '.scores', sep = '\t', index=False)  

if __name__ == "__main__":
    start_time = time.time()
    main()
    print("Finished in "+  str((time.time() - start_time)/60) + " minutes!")

