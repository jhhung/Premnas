import subprocess
import os
import argparse

def Learning_ad_hoc_subpopulation_characteristic():

    print("--------------------------------------------")
    print("Learning ad hoc subpopulation characteristic")
    print("--------------------------------------------")

    # check Premnas images
    all_images = subprocess.check_output('docker images', shell=True).decode("utf-8")
    if all_images.find(PREMNAS_DOCKER) != -1:
        print("Premnas docker exist!")
    else:
        print("Premnas docker doesn't exist")
        print("Download...")
        os.system("docker pull" + PREMNAS_DOCKER + "/" + PREMNAS_DOCKER_VERSION)

    os.system("docker run -ti --rm -v " + args.input_dir + ":/input_dir -v " + args.output_dir + ":/output_dir " + \
            PREMNAS_DOCKER + ":" + PREMNAS_DOCKER_VERSION + " /Premnas.sh " + args.single_cell + " " + args.sc_source)


def Performing_digital_cytometry():

    print("----------------------------")
    print("Performing digital cytometry")
    print("----------------------------")
    
    # check CIBERSOERx images
    all_images = subprocess.check_output('docker images', shell=True).decode("utf-8")
    if all_images.find("cibersortx/fractions") != -1:
        print("CIBERSORTx docker exist!")
    else:
        print("CIBERSORTx docker doesn't exist")
        print("Download...")
        os.system("docker pull cibersortx/fraction") 
    
    if args.mode =! "All":
    	os.system("cp " + args.output_dir + "/subpopulation-characteristic.txt " + args.input_dir)
    	os.system("docker run -v " + args.input_dir + ":/src/data -v " + args.output_dir + ":/src/outdir cibersortx/fractions --username " + \
    		    args.user_name + " --token " + args.token + "--single_cell TRUE --fraction 0 --rmbatchSmode TRUE --perm 500 --refsample " + \
    		    "subpopulation-characteristic.txt --mixture " + args.mixture)
    else:
    	os.system("docker run -v " + args.input_dir + ":/src/data -v " + args.output_dir + ":/src/outdir cibersortx/fractions --username " + \
    		    args.user_name + " --token " + args.token + "--single_cell TRUE --fraction 0 --rmbatchSmode TRUE --perm 500 --refsample " + \
    		    args.cel_characteristic + " --mixture " + args.mixture)



def Analyzing_subpopulation_change():

    print("------------------------------")
    print("Analyzing subpopulation change")
    print("------------------------------")

    os.system("python3 src/Treatment-selection.py " + args.input_dir + "/" + args.metadata + " " + \
            args.output_dir + "/CIBERSORTx_Adjusted.txt -0.9 -0.8") 

    print("Selected cocktail therapy is stored in Treatment-selection-output.csv")

def check_arg(mode):

    if mode != 'Download_CMap_data':
        if args.input_dir is None:
            print("ERROR: input directory should be set")
            return False
        elif args.output_dir is None:
            print("ERROR: output directory should be set")
            return False
    else:
        if args.celltype is None:
            print("ERROR: cell line argument should be set")
            return True

    if mode == "Learn_char" or mode == "All":
        
        if args.single_cell is None:
            print("ERROR: Single cell data must be included")
            return False
        elif args.sc_source is None:
            print("ERROR: Source clone labels must be included for batch correction!")
            return False
    
    elif mode == "Dig_cytomery" or mode == "All":

        if args.user_name is None:
            print("ERROR: Registered user name of COBERSORTx should be set to run the CIBERSORTx images")
            return False
        elif args.token is None:
            print("ERROR: Token of COBERSORTx should be set to run the CIBERSORTx images")
            return False
        elif args.mixture is None:
            print("ERROR: CMap data should be set to run the CIBERSORTx decompostion")
            return False

        if mode == "Dig_cytometry":
            if args.sub_characteristic is None:
                print("ERROR: Single cell characteristic file should be included.")
                return False
    
    elif mode == "Analyze_sub" or mode == "All":

        if args.metadata is None:
            print("ERROR: Information of each perturbation should be set, use Download_CMap_data mode if needed")
            return False

    return True
        

if __name__ == '__main__':

    PREMNAS_DOCKER = "toeric/premnas"
    PREMNAS_DOCKER_VERSION = "1.0.6"
    
    parser = argparse.ArgumentParser(description = "Run Premnas")    
    parser.add_argument('mode', choices=['Learn_char', 'Dig_cytomery', 'Analyze_sub', 'All', 'Download_CMap_data'], help='three main steps of Premnas')
    
    parser.add_argument("-input_dir", help = "Absolute folder path of all your input files")
    parser.add_argument("-output_dir", help = "Absolute folder path of all your output files")

    parser.add_argument("-single_cell", help = "Single cell GEPs data used to learn the subpopulation characteristic")
    parser.add_argument("-sc_source", help = "Source clone labels for all cells")
    
    parser.add_argument("-user_name", help = "Registered user name of CIBERSORTx")
    parser.add_argument("-token", help = "CIBERSORTx container token")
    parser.add_argument("-mixture", help = "CMap data")
    parser.add_argument("-celltype", help = "Cell line name, Options: A375|A549|HCC515|HEPG2|MCF7|PC3|VCAP|HT29")
    parser.add_argument("-sub_characteristic", help = "Subpopulation characteristic file")
    
    parser.add_argument("-metadata", help = "Information of each perturbation experiment in the LINCS L1000 CMap database")
    
    args = parser.parse_args()

    if check_arg(args.mode) == False:
        print("ERROR: Missing important argument!")
        quit()

    if args.mode == 'Learn_char':
        Learning_ad_hoc_subpopulation_characteristic()
    elif args.mode == 'Dig_cytomery':
        Performing_digital_cytometry()
    elif args.mode == 'Analyze_sub':
        Analyzing_subpopulation_change()
    elif args.mode == 'Download_CMap_data':
        os.system("python3 src/generate_L1000_GEP.py -c " + args.celltype)
    elif args.mode == "All":
        Learning_ad_hoc_subpopulation_characteristic()
        Performing_digital_cytometry()
        Analyzing_subpopulation_change()
        print("All DONE!")
    else:
        print("ERRROR: Mode argument need to be Download_CMap_data|Learn_char|Dig_cytomery|Analyze_sub|All")
