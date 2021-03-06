import os
import argparse
import pandas as pd
import numpy as np
import pybedtools
#from .SNPUtility import SNPbedtools
#from .SNPUtility import SNPfeatures
#from .CellTypeAnalyzer import CellTypeAnalysis
#from .SNPMotifAnalyzer import SNPMotifPlot
#from SNPUtility import SNPbedtools
#from SNPUtility import SNPfeatures
#from CellTypeAnalyzer import CellTypeAnalysis
#from SNPMotifAnalyzer import SNPMotifPlot
#from SLDPAnalyzer import RunSLDP

if __package__:
    from .SNPUtility import SNPbedtools
    from .SNPUtility import SNPfeatures
    from .SNPUtility import Preprocess
    from .CellTypeAnalyzer import CellTypeAnalysis
    from .SNPMotifAnalyzer import SNPMotifPlot
    from .SNPScanner import Motif_Scan
else:
    from SNPUtility import SNPbedtools
    from SNPUtility import SNPfeatures
    from SNPUtility import Preprocess
    from CellTypeAnalyzer import CellTypeAnalysis
    from SNPMotifAnalyzer import SNPMotifPlot
    from SNPScanner import Motif_Scan

__version__ = '0.3.0'

def main():
    parser = argparse.ArgumentParser(prog='MotifRaptor', description='Analyze motifs and SNPs in the dataset.')
    subparsers = parser.add_subparsers(help='help for subcommand: celltype, snpmotif, snpfeature, motiffilter, motifspecific, snpspecific', dest="command")
    parser.add_argument('--version', action='version', version='%(prog)s 0.3.0')

    parser_k = subparsers.add_parser('preprocess', help='Pre-process the summary statistics')
    parser_k.add_argument('-gwas', '--gwas_summary', type=str, help='GWAS summary statistic file',dest="summarystatsfile")
    parser_k.add_argument('-cn', '--column_numbers', type=str, help='provide six column numbers for information in such order: ID,CHR,POS,A1,A2,SCORE eg. 1,2,3,4,5,6',dest="column_numbers")
    parser_k.add_argument('-st', '--score_type', type=str, default='pvalue',help='Score type in GWAS summary statistic file: pvalue or zscore or chisquare',dest="score_type")
    parser_k.add_argument('-th', '--threshold', type=float, default=5E-8,help='threads for pvalue - default 5E-8',dest="score_pvalue_threshold")

    parser_l = subparsers.add_parser('preprocess_ukbb_v3', help='Pre-process the summary statistics from UKBB version 3 TSV files')
    parser_l.add_argument('-gwas', '--gwas_summary', type=str, help='GWAS summary statistic - UKBB v3 tsv filename',dest="summarystatsfile")
    parser_l.add_argument('-th', '--threshold', type=float, default=5E-8,help='threads for pvalue - default 5E-8',dest="score_pvalue_threshold")

    parser_a = subparsers.add_parser('celltype', help='cell type or tissue type analysis help')
    parser_a.add_argument('-vcf', '--snp_hit_withseq', type=str, help='need header and columns in this text file with sequence (CHR is only a number): ID	CHR	POS	REF	ALT',dest="snp_hit_vcf")
    parser_a.add_argument('-sh', '--snp_hit', type=str, help='need header and columns in this text file (CHR is only a number): ID CHR     POS',dest="snp_hit")
    parser_a.add_argument('-snh', '--snp_non_hit', type=str, help='need header and columns in this text file (CHR is only a number): ID CHR     POS',dest="snp_non_hit")    
    parser_a.add_argument('-wd', '--workdir', type=str, help='Working directory',dest="workdir")
    parser_a.add_argument('-p', '--threads', type=int, help='number of threads',dest="thread_num")

    parser_b = subparsers.add_parser('snpmotif', help='snp motif test help')
    parser_b.add_argument('-wd', '--workdir', type=str, help='Working directory',dest="workdir")
    parser_b.add_argument('-c', '--cell_type', type=str, help='Cell type or Tissue type ID',dest="cell_type")
    parser_b.add_argument('-sb', '--snp_hit_bed', type=str, help='hit snps on union bed file, usually from step1',dest="hit_snp_union_bed")
    parser_b.add_argument('-se', '--snp_hit_seq', type=str, help='hit snps with sequence information, usually from step1',dest="hit_snp_union")
    parser_b.add_argument('-bg', '--bg_snp', type=str, help='background snp list file or (genome)',dest="bg_snps")
    parser_b.add_argument('-m', '--motifs', type=str, help='motif list file, no header, or (all)',dest="motif_list")
    parser_b.add_argument('-p', '--threads', type=int, help='number of threads',dest="thread_num")

    parser_m = subparsers.add_parser('gkm_snpmotif', help='gkm snp motif test help')
    parser_m.add_argument('-wd', '--workdir', type=str, help='Working directory',dest="workdir")
    parser_m.add_argument('-c', '--cell_type', type=str, help='Cell type or Tissue type ID',dest="cell_type")
    parser_m.add_argument('-sb', '--snp_hit_bed', type=str, help='hit snps on union bed file, usually from step1',dest="hit_snp_union_bed")
    parser_m.add_argument('-se', '--snp_hit_seq', type=str, help='hit snps with sequence information, usually from step1',dest="hit_snp_union")
    parser_m.add_argument('-bg', '--bg_snp', type=str, help='background snp list file or (genome)',dest="bg_snps")
    parser_m.add_argument('-m', '--motifs', type=str, help='motif list file, no header, or (all)',dest="motif_list")
    parser_m.add_argument('-p', '--threads', type=int, help='number of threads',dest="thread_num")

    parser_d = subparsers.add_parser('snpfeature', help='snp feature help')
    parser_d.add_argument('-wd', '--workdir', type=str, help='Working directory',dest="workdir")
    parser_d.add_argument('-c', '--cell_type', type=str, help='Cell type or Tissue type ID',dest="cell_type")
    parser_d.add_argument('-cb', '--snp_bed_files', type=str, help='SNP cell type bed file folder, usually from step1',dest="snp_bed_files")


    parser_c = subparsers.add_parser('motiffilter', help='motifs filtering help')
    parser_c.add_argument('-wd', '--workdir', type=str, help='Working directory',dest="workdir")
    parser_c.add_argument('-ms', '--motif_summary', type=str, help='Motif Summary File, usually from step2',dest="motiffile")
    
    parser_d = subparsers.add_parser('motifspecific', help='motifs specific analysis help')
    parser_d.add_argument('-wd', '--workdir', type=str, help='Working directory',dest="workdir")
    parser_d.add_argument('-sm', '--snp_motif_file', type=str, help='SNP motif pair-wise list File, usually from step2',dest="snp_motif_file")
    parser_d.add_argument('-md', '--motif_id', type=str, help='motif id with name, in the format of motifID__NAME',dest="motif_id_name")
    parser_d.add_argument('-bs', '--bg_score_folder', type=str, help='background score folder, usually from step2',dest="bg_score_folder")

    parser_e = subparsers.add_parser('snpspecific', help='SNP specific analysis help')
    parser_e.add_argument('-wd', '--workdir', type=str, help='Working directory',dest="workdir")
    parser_e.add_argument('-sm', '--snp_motif_file', type=str, help='SNP motif pair-wise list File, usually from step2',dest="snp_motif_file")
    parser_e.add_argument('-snp', '--snp_id', type=str, help='SNP id',dest="snp_id")

    parser_f = subparsers.add_parser('snpmotifradar', help='SNP motif radar plot help')
    parser_f.add_argument('-wd', '--workdir', type=str, help='Working directory',dest="workdir")
    parser_f.add_argument('-sm', '--snp_motif_file', type=str, help='SNP motif pair-wise list File, usually from step2',dest="snp_motif_file")
    parser_f.add_argument('-sf', '--snp_feature_file', type=str, help='SNP feature file, usually from step2',dest="snp_feature_file")
    parser_f.add_argument('-pid', '--snp_motif_id', type=str, help='SNP motif pair-wise ID',dest="snp_motif_id")


    parser_j = subparsers.add_parser('snpindex', help='index the SNPs (with flanking sequences) help')
    parser_j.add_argument('-vcf', '--vcf_filename', type=str, help='input - VCF file for SNPs, first five columns need to be CHR,POS,ID,REF,ALT',dest="vcf_filename")
    parser_j.add_argument('-gi', '--indexed_genome_db', type=str, help='output - indexed genome_database_folder',dest="indexed_genome_db")
    parser_j.add_argument('-mk', '--mksary', type=str, action='store', default='', help='Mksary program path',dest="mksary_path")
    parser_j.add_argument('-p', '--threads', type=int, help='number of threads',dest="thread_num")

    parser_g = subparsers.add_parser('snpscan', help='scan SNP database (already indexed) help')
    #parser_g.add_argument('-g', '--genome_db', type=str, help='genome_database_folder',dest="genome_db")
    parser_g.add_argument('-gi', '--indexed_genome_db', type=str, help='indexed genome_database_folder',dest="indexed_genome_db")
    parser_g.add_argument('-pfm', '--pfm_folder', type=str, help='motif pmf files folder',dest="pfm_folder")
    parser_g.add_argument('-mo', '--motifscan_output', type=str, action='store', default='./motifscanfiles', help='Motif Scan Ouput Folder',dest="outputmotifscandir")
    parser_g.add_argument('-p', '--threads', type=int, help='number of threads',dest="thread_num")

    parser_h = subparsers.add_parser('gkm_snpscan', help='scan SNP database (using gkmSVM models) help')
    parser_h.add_argument('-b', '--background', type=str, help='background snp filename',dest="background_snp")    
    parser_h.add_argument('-m', '--models', type=str, help='models names',dest="models")    
    parser_h.add_argument('-o', '--motifscan_output', type=str, action='store', default='./motifscanfiles', help='Motif Scan Ouput Folder',dest="outputmotifscandir")
    parser_h.add_argument('-p', '--threads', type=int, help='number of threads',dest="thread_num")

    parser_i = subparsers.add_parser('set', help='Set Path and Global Values')
    parser_i.add_argument('-pn', '--parametername', type=str, help='Parameter Name',dest="parametername")
    parser_i.add_argument('-pv', '--parametervalue', type=str, help='Parameter Value',dest="parametervalue")

    parser_l = subparsers.add_parser('info', help='Get Informationa and Print Global Values')

    package_path=os.path.dirname(os.path.abspath(__file__))
    globalfile=os.path.join(package_path,"globalpara.txt")
    p_path_1=os.path.join(package_path,"Database/hg19/")
    p_path_2=os.path.join(package_path,"Database/hg19/motifdatabase/") 
    global_df=pd.DataFrame({'para':['databasedir','motifdatabasedir'],
                            'value':[p_path_1,p_path_2]})
    if os.path.exists(globalfile):
        global_df=pd.read_csv(globalfile,header=0,sep="\t")
    global_df.to_csv(globalfile,sep="\t",header=True,index=None)

    args = parser.parse_args()
    print(args) 
    #parser.print_help()
    #package_path=os.path.dirname(__file__)
    global_df=pd.read_csv(globalfile,header=0,sep="\t")
    package_path_1 = global_df.loc[global_df['para']=="databasedir"]['value'].values[0]
    package_path_2 = global_df.loc[global_df['para']=="motifdatabasedir"]['value'].values[0]

    if args.command=="celltype":
        print("Command: fetch the snp information ...")
        if not os.path.exists(args.workdir):
            os.mkdir(args.workdir)
        genomefilename=os.path.join(package_path_1,"hg19.2bit")
        vcftoseq_command="python3 "+os.path.join(package_path,"SNPUtility/VCFtoSeq.py")+" "+args.snp_hit_vcf+" "+genomefilename+" 0,1,2,3,4 ID,CHR,POS,REF,ALT 30 "+str(args.thread_num)+" "+os.path.join(args.workdir,"hit_snps_seq_pickle")
        print("Command: "+vcftoseq_command)
        status = os.system(vcftoseq_command)
        print("Command: merge with union DHS")
        DHS_union_filename=os.path.join(package_path_1,"union_DHS.bed")  
        ####real SNP hits need to calculate on the fly in later steps      
        SNPbedtools.generateSNPBed_fast(os.path.join(args.workdir,"hit_snps_seq_pickle.df.txt"),os.path.join(args.workdir,"tempSNP_hits_sorted.bed"))
        SNPbedtools.overlapSNPandBed(os.path.join(args.workdir,"tempSNP_hits_sorted.bed"), DHS_union_filename, os.path.join(args.workdir,"hitSNP_DHSunion_list.bed"))
        ####only for this cell type purpose
        SNPbedtools.generateSNPBed_fast(args.snp_hit,os.path.join(args.workdir,"allSNP_hits_sorted.bed"))
        SNPbedtools.overlapSNPandBed(os.path.join(args.workdir,"allSNP_hits_sorted.bed"), DHS_union_filename, os.path.join(args.workdir,"allhitSNP_DHSunion.bed"))
        SNPbedtools.generateSNPBed_fast(args.snp_non_hit,os.path.join(args.workdir,"allSNP_nonhits_sorted.bed"))
        SNPbedtools.overlapSNPandBed(os.path.join(args.workdir,"allSNP_nonhits_sorted.bed"), DHS_union_filename, os.path.join(args.workdir,"allnonhitSNP_DHSunion.bed"))
        print("Command: analyze cell type or tissue type ...")
        hit_SNP_df=pd.read_csv(args.snp_hit,sep='\t',header=0)
        hit_SNP_df=hit_SNP_df.set_index('ID')
        non_hit_SNP_df=pd.read_csv(args.snp_non_hit,sep='\t',header=0)
        non_hit_SNP_df=non_hit_SNP_df.set_index('ID')
        CellTypeAnalysis.celltypeanalysis_main(package_path_1,args.workdir,hit_SNP_df, non_hit_SNP_df, args.thread_num)
        print("Command: plot figures ...")
        print("Done!!!")

    elif args.command=="snpmotif":
        if not os.path.exists(args.workdir):
            os.mkdir(args.workdir)
        celltype=args.cell_type.replace("%20"," ")
        hit_snp_union_bed=args.hit_snp_union_bed
        hit_snp_union_file=args.hit_snp_union
        bg_genome=args.bg_snps
        print("Command: prepare cell type information...")
        celltypemappingfile=os.path.join(package_path_1,"download_meta_tier123.txt")
        celltypetable=pd.read_csv(celltypemappingfile,sep='\t')
        fileaccession=celltypetable[(celltypetable['Biosample term name']==celltype)&(celltypetable['Assay']=='DNase-seq')]['File accession'].iloc[0]
        expressionaccession=celltypetable[(celltypetable['Biosample term name']==celltype)&(celltypetable['Assay']=='RNA-seq')]['File accession'].iloc[0]    
        hit_snp_union = pybedtools.BedTool(hit_snp_union_bed)
        DHSdir=os.path.join(package_path_1,"DHS_build/")
        targetbedfile=os.path.join(DHSdir, fileaccession+'_target.bed')
        print("File accession :"+fileaccession)
        print("Expression : "+expressionaccession) 
        b = pybedtools.BedTool(targetbedfile)
        tissue_id=fileaccession
        (hit_snp_union+b).moveto(os.path.join(args.workdir,"tempSNP_"+tissue_id+".txt"))

        hit_snps_tissue=pd.read_csv(os.path.join(args.workdir,"tempSNP_"+tissue_id+".txt"),header=None,sep="\t")
        hit_snps=pd.read_csv(hit_snp_union_file,header=0,sep="\t")
        hit_snps=hit_snps.set_index('ID')
        hit_snps_sub=hit_snps.loc[hit_snps_tissue[3]]
        hit_snps_sub.to_csv(os.path.join(args.workdir,"hit_snps_seq_pickle_"+tissue_id+".df.txt"), sep='\t')

        print("Command: calculating snp and motif disruptions...")
        target_snp_filename=os.path.join(args.workdir,"hit_snps_seq_pickle_"+tissue_id+".df.txt")#eg.hit_snps_seq_pickle_ENCFF512IML.df.txt
        background_snp_filename=bg_genome
        motif_list_filename=args.motif_list
        motif_expression_folder=os.path.join(package_path_1,"TF_expression")
        motif_expression_filename=os.path.join(motif_expression_folder,expressionaccession+".tsv")
        motif_pfm_folder=os.path.join(package_path_2,"pfmfiles")
        motif_scan_folder=os.path.join(package_path_2,"motifscanfiles")
        output_filename=os.path.join(args.workdir,"result_new_df_motifs_"+tissue_id+".txt")
        num_of_threads=args.thread_num
   
        snpmotifscan_command="python3 "+os.path.join(package_path,"SNPMotifAnalyzer/SNPMotifScan.py pfmscan") \
	    +" -t "+target_snp_filename \
	    +" -bg "+background_snp_filename \
	    +" -m "+motif_list_filename \
	    +" -e "+motif_expression_filename \
	    +" -pfm "+motif_pfm_folder \
	    +" -score "+motif_scan_folder \
	    +" -mo "+output_filename \
	    +" -p "+str(num_of_threads) 
        print("Command : snpmotifscan_command")
        status=os.system(snpmotifscan_command)
        
        print("Command : running test on motifs...")
        snp_motif_filename=output_filename
        snp_background_filename=bg_genome
        #motif_expression_filename="ENCFF557VGS.tsv"
        #motif_scan_folder="motifscanfiles"
        outputdir=os.path.join(args.workdir,"motif_result")

        #num_of_threads=40
        if not os.path.exists(outputdir):
            os.mkdir(outputdir)
        snpmotifanalysis="python3 "+os.path.join(package_path,"SNPMotifAnalyzer/SNPMotifAnalysis.py motif_test") \
	    +" -r "+snp_motif_filename \
	    +" -bg "+snp_background_filename \
	    +" -e "+motif_expression_filename \
	    +" -score "+motif_scan_folder \
	    +" -od "+outputdir \
	    +" -p "+str(num_of_threads)
        print("Command : "+snpmotifanalysis)
        status=os.system(snpmotifanalysis)
        print("Done!!!")

    elif args.command=="gkm_snpmotif":
        if not os.path.exists(args.workdir):
            os.mkdir(args.workdir)
        celltype=args.cell_type.replace("%20"," ")
        hit_snp_union_bed=args.hit_snp_union_bed
        hit_snp_union_file=args.hit_snp_union
        bg_genome=args.bg_snps
        print("Command: prepare cell type information...")
        celltypemappingfile=os.path.join(package_path_1,"download_meta_tier123.txt")
        celltypetable=pd.read_csv(celltypemappingfile,sep='\t')
        fileaccession=celltypetable[(celltypetable['Biosample term name']==celltype)&(celltypetable['Assay']=='DNase-seq')]['File accession'].iloc[0]
        expressionaccession=celltypetable[(celltypetable['Biosample term name']==celltype)&(celltypetable['Assay']=='RNA-seq')]['File accession'].iloc[0]    
        hit_snp_union = pybedtools.BedTool(hit_snp_union_bed)
        DHSdir=os.path.join(package_path_1,"DHS_build/")
        targetbedfile=os.path.join(DHSdir, fileaccession+'_target.bed')
        print("File accession :"+fileaccession)
        print("Expression : "+expressionaccession) 
        b = pybedtools.BedTool(targetbedfile)
        tissue_id=fileaccession
        (hit_snp_union+b).moveto(os.path.join(args.workdir,"tempSNP_"+tissue_id+".txt"))

        hit_snps_tissue=pd.read_csv(os.path.join(args.workdir,"tempSNP_"+tissue_id+".txt"),header=None,sep="\t")
        hit_snps=pd.read_csv(hit_snp_union_file,header=0,sep="\t")
        hit_snps=hit_snps.set_index('ID')
        hit_snps_sub=hit_snps.loc[hit_snps_tissue[3]]
        hit_snps_sub.to_csv(os.path.join(args.workdir,"hit_snps_seq_pickle_"+tissue_id+".df.txt"), sep='\t')

        print("Command: calculating snp and motif disruptions...")
        target_snp_filename=os.path.join(args.workdir,"hit_snps_seq_pickle_"+tissue_id+".df.txt")#eg.hit_snps_seq_pickle_ENCFF512IML.df.txt
        background_snp_filename=bg_genome
        motif_list_filename=args.motif_list
        motif_expression_folder=os.path.join(package_path_1,"TF_expression")
        motif_expression_filename=os.path.join(motif_expression_folder,expressionaccession+".tsv")
        gkm_score_folder=os.path.join(package_path_2,"gkm_score")
        motif_scan_folder=os.path.join(package_path_2,"motifscanfiles")
        output_filename=os.path.join(args.workdir,"result_new_df_motifs_"+tissue_id+".txt")
        num_of_threads=args.thread_num

        snpmotifscan_command="python3 "+os.path.join(package_path,"SNPMotifAnalyzer/SNPMotifScan.py gkm_scan") \
        +" -t "+target_snp_filename \
        +" -bg "+background_snp_filename \
        +" -m "+motif_list_filename \
        +" -e "+motif_expression_filename \
        +" -gkm "+gkm_score_folder \
        +" -score "+motif_scan_folder \
        +" -mo "+output_filename \
        +" -p "+str(num_of_threads) 
        print("Command : snpmotifscan_command")
        status=os.system(snpmotifscan_command)
        
        print("Command : running test on motifs...")
        snp_motif_filename=output_filename
        snp_background_filename=bg_genome
        #motif_expression_filename="ENCFF557VGS.tsv"
        #motif_scan_folder="motifscanfiles"
        outputdir=os.path.join(args.workdir,"motif_result")

        #num_of_threads=40
        if not os.path.exists(outputdir):
            os.mkdir(outputdir)
        snpmotifanalysis="python3 "+os.path.join(package_path,"SNPMotifAnalyzer/SNPMotifAnalysis.py motif_test") \
        +" -r "+snp_motif_filename \
        +" -bg "+snp_background_filename \
        +" -e "+motif_expression_filename \
        +" -score "+motif_scan_folder \
        +" -od "+outputdir \
        +" -p "+str(num_of_threads)
        print("Command : "+snpmotifanalysis)
        status=os.system(snpmotifanalysis)
        print("Done!!!")

    elif args.command=="snpfeature":
        bedfiledir=args.snp_bed_files
        if not os.path.exists(args.workdir):
            os.mkdir(args.workdir)
        workdir=args.workdir
        celltype2=args.cell_type.replace("%20"," ")
        conservationfolder=os.path.join(package_path_1,"SNP_conservation")
        bed_hit=os.path.join(bedfiledir,celltype2+".result.bed")
        bed_nonhit=os.path.join(bedfiledir,celltype2+".result_bg.bed")
        catoutput=os.path.join(args.workdir,"tempallSNP.bed")
        catcommand="cat \""+bed_hit+"\" \""+bed_nonhit+"\">"+catoutput
        os.system(catcommand)
        a = pybedtools.BedTool(catoutput)
        celltype=args.cell_type.replace("%20"," ")
        celltypemappingfile=os.path.join(package_path_1,"download_meta_tier123.txt")
        celltypetable=pd.read_csv(celltypemappingfile,sep='\t')
        fileaccession=celltypetable[(celltypetable['Biosample term name']==celltype)&(celltypetable['Assay']=='DNase-seq')]['File accession'].iloc[0]
        DHSdir=os.path.join(package_path_1,"DHS_build/")
        targetbedfile=os.path.join(DHSdir, fileaccession+'_target.bed')
        b = pybedtools.BedTool(targetbedfile)
        tissue_id=fileaccession
        (a+b).moveto(os.path.join(workdir,"tempallSNP_"+tissue_id+".bed"))

        #remove exons
        snpbed=os.path.join(workdir,"tempallSNP_"+tissue_id+".bed")
        exonbed=os.path.join(package_path_1,"Homo_sapiens_exons.GRCh37.75.bed")
        a = pybedtools.BedTool(snpbed)
        a=a.sort()
        b=pybedtools.BedTool(exonbed)
        b=b.sort()
        c = a.intersect(b, v=True)
        snpnoexonbed=os.path.join(workdir,"tempallSNP_"+tissue_id+".noexon.bed")
        c.moveto(snpnoexonbed)

        featurefolder=conservationfolder
        SNPbedfile=snpbed
        SNPfeaturedf_outfile=os.path.join(workdir,"SNP_"+tissue_id+"_features.txt")
        snpnoexon_df=pd.read_csv(snpnoexonbed, sep="\t", header=None)
        SNPnormalize_list=snpnoexon_df[3]
        SNPfeatures.SNPfeatures_main(featurefolder, SNPbedfile, SNPfeaturedf_outfile, SNPnormalize_list)
        print("Done!!!")

    elif args.command=="motiffilter":
        motiffile=args.motiffile
        all_expression_file=os.path.join(package_path_1,"Expression_all/all_expression.csv")
        if not os.path.exists(args.workdir):
            os.mkdir(args.workdir)
        motiffigure_filename=os.path.join(args.workdir,"plot_motif_summary.pdf")
        output_filename=os.path.join(args.workdir,"motif_summary.txt")
        SNPMotifPlot.plot_motif_scattering_main(motiffile,all_expression_file,output_filename,motiffigure_filename,0.8, 1,0,0.1)
        print("Done!!!")
    
    elif args.command=="motifspecific":
        if not os.path.exists(args.workdir):
            os.mkdir(args.workdir)
        snp_motif_result_file=args.snp_motif_file
        bg_score_folder=args.bg_score_folder
        motif_id_name=args.motif_id_name
        pngfilename=os.path.join(args.workdir,"plot_distribution_"+motif_id_name+".png")
        SNPMotifPlot.plot_motif_specific_main(snp_motif_result_file, bg_score_folder, motif_id_name, pngfilename)
        print("Done!!!")

    elif args.command=="snpspecific":
        if not os.path.exists(args.workdir):
            os.mkdir(args.workdir)
        snp_motif_result_file=args.snp_motif_file
        snp_id=args.snp_id
        pdffilename=os.path.join(args.workdir,"plot_TF_for_"+snp_id+".pdf")
        all_expression_file=os.path.join(package_path_1,"Expression_all/all_expression.csv")
        SNPMotifPlot.plot_snp_specific_main(snp_motif_result_file, snp_id, all_expression_file,pdffilename)
        print("Done!!!")

    elif args.command=="snpmotifradar":
        if not os.path.exists(args.workdir):
            os.mkdir(args.workdir)
        snp_motif_result_file=args.snp_motif_file
        rsid_motifid=args.snp_motif_id
        rsid_motifid_list=rsid_motifid.split(":")
        snp_id=rsid_motifid_list[0]
        motif_id_name=rsid_motifid_list[1]
        print("SNP : "+snp_id)
        print("motif :"+motif_id_name)
        pdffilename=os.path.join(args.workdir,"plot_radar_for_"+snp_id+"_"+motif_id_name+".pdf")
        all_expression_file=os.path.join(package_path_1,"Expression_all/all_expression.csv")
        #conservation_folder=os.path.join(package_path,"Database/hg19/SNP_conservation")
        #SNPMotifPlot.plot_motif_snp_pair_main(snp_motif_result_file, all_expression_file,snp_id, motif_id_name, rsid_motifid, pdffilename, conservation_folder)
        catofile=os.path.join(package_path_1,"SNP_Motif_score/CATO.txt")
        snpfeaturefile=args.snp_feature_file
        SNPMotifPlot.plot_motif_snp_pair_main(snp_motif_result_file, all_expression_file,snp_id, motif_id_name, rsid_motifid, pdffilename, snpfeaturefile,catofile)
        print("Done!!!")

    elif args.command=="preprocess":
        sumstatsfile = args.summarystatsfile
        score_type = args.score_type
        pthreshold = args.score_pvalue_threshold

        usecolumnsstring=args.column_numbers.strip().split(',')
        usecolumns_list=[int(ai)-1 for ai in usecolumnsstring]
        Preprocess.preprocess_summary_statistics_general(sumstatsfile,score_type,pthreshold,usecolumns_list)

        print("Done!!!")

    elif args.command=="preprocess_ukbb_v3":
        ukbb_tsv_filename = args.summarystatsfile
        
        pthreshold = args.score_pvalue_threshold

        Preprocess.preprocess_from_ukbb_v3(ukbb_tsv_filename,pthreshold)
        print("Done!!!")

    elif args.command=="snpindex":
        outfolder=args.indexed_genome_db
        if not os.path.exists(outfolder):
            os.mkdir(outfolder)
        vcf_filename = args.vcf_filename
        usecolumns_list=[0,1,2,3,4]
        column_renames_list=['CHR','POS','ID','REF','ALT']
        window_size=30
        genomefilename=os.path.join(package_path_1,"hg19.2bit")
        Preprocess.preprocess_vcf_for_scan(genomefilename, vcf_filename, usecolumns_list, column_renames_list, outfolder, window_size) 

        mksary_path=args.mksary_path
        if mksary_path=="":
            mksary_path = os.path.join(package_path,"SNPScanner")
        genome_db = outfolder
        Motif_Scan.run_mksary(mksary_path, genome_db, genome_db, args.thread_num)
        print("Done!!!")

    elif args.command=="snpscan":
        indexed_genome_db=args.indexed_genome_db
        pfm_folder=args.pfm_folder
        outputmotifscandir=args.outputmotifscandir
        Motif_Scan.run_scan_motif(indexed_genome_db, pfm_folder, outputmotifscandir, args.thread_num)

        snpmotifscale_command="python3 "+os.path.join(package_path,"SNPMotifAnalyzer/SNPMotifScan.py scale") \
            +" -pfm "+pfm_folder \
            +" -score "+outputmotifscandir \
            +" -p "+str(args.thread_num)
        print("Command : snpmotifscale_command")
        status=os.system(snpmotifscale_command)

        print("Done!!!")

    elif args.command=="gkm_snpscan":
        background_snp_filename=args.background_snp
        gkm_score_folder=os.path.join(package_path_2,"gkm_score")
        models=args.models
        outputmotifscandir=args.outputmotifscandir
        
        snpmotifscore_command="python3 "+os.path.join(package_path,"SNPMotifAnalyzer/SNPMotifScore.py") \
        +" -b "+background_snp_filename \
        +" -wd "+gkm_score_folder \
        +" -m "+models \
        +" -o "+ outputmotifscandir
        print("Command: snpmotifscore_command")
        status=os.system(snpmotifscore_command)
        print("Done!!!")


        snpmotifscale_command="python3 "+os.path.join(package_path,"SNPMotifAnalyzer/SNPMotifScan.py gkm_scale") \
            +" -s " + outputmotifscandir\
            +" -mo " + outputmotifscandir\
            +" -p "+str(args.thread_num)
        print("Command: snpmotifscale_command")
        status=os.system(snpmotifscale_command)

        print("Done!!!")

    elif args.command=="set":
        parametername=args.parametername
        parametervalue=args.parametervalue
        if parametername in global_df['para'].tolist():
            global_df.loc[global_df['para'] == parametername, 'value'] = parametervalue
            print("Set "+parametername+ " As "+parametervalue)
        global_df.to_csv(globalfile,sep="\t",header=True,index=None)
        print("Done!!!")
    elif args.command=="info":
        print("Motif-Raptor version: "+__version__)
        print("package path: "+package_path)
        print("Global Parameters : ")
        print(global_df)
if __name__=="__main__":
    main()
