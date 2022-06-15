import os
import sys
import re
import glob
import shutil
import subprocess

print("Program to analyse:\nGaussian = 1\nPsi4 = 2\nOrca = 3\nCrest = 4\nGAMESS = 5\nMisc Functions = 6")
program_type = str(input())

class CrestClass:
    def __init__(self):
        #unpackfile = "crest_conformers.xyz"
        #enercutoff = 3
        #molname = "test"
        print("File to unpack:")
        unpackfile = input()
        print("Energy Cutoff:")
        enercutoff = input()        
        print("Molecule Name:")
        molname = input()
        self.unpackfile = unpackfile
        self.enercutoff = enercutoff        
        self.molname = molname
    def unpack(self):
        with open(self.unpackfile, "r") as f:
            allcoords = f.read()
            coordlines = allcoords.split("\n")
            atom_numbers = coordlines[0].strip()
            seperated_coords = re.sub("\s{2,}" + atom_numbers + "\s*\n","\nnewmolecule\n",allcoords).split("newmolecule")
            del seperated_coords[0]
        if not os.path.exists("conformers"):
            os.mkdir("conformers")
            os.chdir("conformers")
        else:
            os.chdir("conformers")
            cwd = os.getcwd()
            for f in os.listdir(cwd):
                os.remove(os.path.join(cwd, f))
        base_energy = 0
        for index, coord in enumerate(seperated_coords):
            if index == 0:
                base_energy = float(coord.split()[0])
            energy_difference = (float(coord.split()[0]) - base_energy) * 2625.5
            if  energy_difference <= float(self.enercutoff):
                with open(f'{self.molname}-{str(index)}.xyz', "w") as f:
                    f.write("  " + atom_numbers)
                    f.write(coord)
                    print(f"Saving file {self.molname}-{str(index)}|{energy_difference}")
            else:
                pass

        return

class Psi4Class:
    def __init__(self):
        cwd = os.getcwd()
        self.cwd = cwd
    def createjob(self):
        with open("job.template", "r") as f:
            job_template = f.read()
        with open("input.template", "r") as f:
            input_template = f.read()
            template_frags = re.split("[a-zA-Z]{1,2}\s{5,}-{0,1}\d{1,3}.\d*\s{5,}-{0,1}\d{1,3}.\d*\s{5,}-{0,1}\d{1,3}.\d*", input_template)
            jobfilepart1 = template_frags[0]
            jobfilepart2 = template_frags[-1]
            xyzlist = glob.glob("*xyz")
            for xyz in xyzlist:
                molname = xyz.split(".xyz")[0]
                job_file = job_template.replace("job.template",molname)
                with open(xyz, "r") as f:
                    xyz_coords = f.readlines()[2:]
                if not os.path.exists(molname):       
                    os.mkdir(molname)
                    os.chdir(molname)
                else:
                    os.chdir(molname)                 
                with open(molname + ".job", "w") as f:
                    f.write(job_file)
                with open(molname + ".inp", "w") as f:
                    f.write(jobfilepart1+ "".join(xyz_coords) + jobfilepart2)
                os.chdir(self.cwd)
                shutil.copyfile(xyz, f'{molname}/{xyz}')
    def calcenergy(self):
        log_file_list = []
        energies_list = []
        for dirpath, dirnames, filenames in os.walk("."):
            for filename in [f for f in filenames if f.endswith(".out")]:
                log_file_list.append(os.path.join(dirpath, filename))
        if not "psi4energies.csv" in glob.glob("*.csv"):
            pass
        else:
            os.remove("psi4energies.csv")
        for log_file in log_file_list:
            totalenergy = "Not Found"
            refenergy = "Not Found"
            ssenergy = "Not Found"
            osenergy = "Not Found"
            if "ccpvdz" in log_file:
                basis_set = "cc-pVDZ"
            elif "ccpvtz" in log_file:
                basis_set = "cc-pVTZ"
            elif "ccpvqz" in log_file:
                basis_set = "cc-pVQZ"
            else:
                basis_set = "Not Found"
            with open(log_file, "r") as f:
                log_lines = f.readlines()
                for line in log_lines:
                    if "Total Energy" in line and "[Eh]" in line and not "SCS" in line:
                        totalenergy = line.split()[3]
                    if "Reference Energy" in line and not "SCS" in line:
                        refenergy = line.split()[3]
                    if "Same-Spin Energy" in line and not "SCS" in line:
                        ssenergy = line.split()[3]
                    if "Opposite-Spin Energy" in line and not "SCS" in line:
                        osenergy = line.split()[3]
                correnergy = ssenergy + osenergy
                energies_list.append(f'{log_file},{basis_set},{totalenergy},{refenergy},{ssenergy},{osenergy},{correnergy}')
        with open("psi4energies.csv", "w") as f:
            f.write("Molecule,Basis Set,Total Energy,Reference Energy,SS Energy,OS Energy,Correlation Energy\n")
            for energy in energies_list:
                f.write(f'{energy}\n')
class GaussianClass:
    def createjob(self):
        with open("job.template", "r") as f:
            template = f.read()
            template_frags = re.split("[a-zA-Z]{1,2}\s*-{0,1}\d{1,3}\.\d*\s*-{0,1}\d{1,3}\.\d*\s*-{0,1}\d{1,3}\.\d*", template)
            jobfilepart1 = template_frags[0]
            jobfilepart2 = template_frags[-1]
        xyzlist = glob.glob("*xyz")
        for xyz in xyzlist:
            molname = xyz.split(".xyz")[0]
            molname = molname.replace('_convError', '')
            print(molname)
            jobfilepart1 = jobfilepart1.replace("job.template",f'{molname}')
            with open(xyz, "r") as f:
                xyz_coords = f.readlines()[2:]
            with open(molname + ".job", "w") as f:
                print(molname)
                f.write(jobfilepart1+ "".join(xyz_coords) + jobfilepart2)
            jobfilepart1 = jobfilepart1.replace(f'{molname}',"job.template")
    def processjob(self):
        job_file_list = glob.glob('*.job')
        cwd = os.getcwd()
        job_outputfile_list = glob.glob("*.job.o*")
        qstat = subprocess.check_output(["qstat"])
        atom_number_dict = {
            "1": "H",
            "6": "C",
            "7": "N",
            "5": "B",
            "8": "O",
            "9": "F",
            "15": "P",
            "16": "S",
            "17": "Cl"
            }
        for job_file in job_file_list:
            job_file_segments = job_file.split(".")
            molecule_name = job_file_segments[0]
        qstat = qstat.decode("utf=8")
        lines = qstat.splitlines()[2:]
        job_numbers_submitted = []
        p = False
        for l in lines:
            job = l.split()
            job_number = job[0]
            job_numbers_submitted.append(job_number)
        for job_outputfile in job_outputfile_list:
            job_outputfile_segments = job_outputfile.split(".o")
            job_outputfile_name_segments = job_outputfile_segments[0].split(".")
            molecule_output_name = job_outputfile_name_segments[0]
            job_ID = job_outputfile_segments[1]
            if job_ID in job_numbers_submitted:
                job_allocation = "still going"
            else:
                job_allocation = "finished"
            log_file = molecule_output_name + ".log"
            chk_file = molecule_output_name + ".chk"
            po_file = molecule_output_name + ".job.po" + job_ID
            job_file = molecule_output_name + ".job"
            print(log_file)
            if "finished" in job_allocation:
                xyz_coords = [] # save coordinates into this list
                xyz_original = []
                xyz_reformatted = []
                with open(log_file) as f:
                    ff2 = f.read()
                    f.seek(0)
                    ff = f.readlines()
                    if "Normal termination" in ff2:
                        job_allocation = "Completed Successfully"
                    #elif "Erroneous read" in ff2:
                    #    job_allocation = "Read Error"
                    elif "Error termination request processed by link 9999." in ff2:
                        job_allocation = "Maximum iterations reached"
                    elif "galloc:  could not allocate memory." in ff2:
                        job_allocation = "Memory allocation failure"
                    else:
                        job_allocation = "Unknown error"
                    rf = list(reversed(ff))
                    if "Read Error" in job_allocation:
                        for l in rf:
                            if "Rotational constants (GHZ)" in l or "Cartesian Forces" in l:
                                p = True # results from last iteration begin
                            elif "Number     Number       Type" in l:
                                p = False # results from last iteration end
                                break # don't bother reading further
                            if p:
                                xyz_coords.append(l.strip()) # save results
                        xyz_coords = list(reversed(xyz_coords))
                        xyz_coords = xyz_coords[1:-2] #cutting results
                        for xyz in xyz_coords:
                            xyz = xyz.split()
                            xyz[1] = atom_number_dict.get(xyz[1])
                            del xyz[0]
                            del xyz[1]
                            xyz = "   ".join(xyz)
                            xyz_reformatted.append(xyz)
                        with open(job_file,"r+") as f:
                            ff = f.readlines() #needed for replacing later on
                            f.seek(0) # to get f back to start
                            for line in f:
                                if re.match("^\s*\d{1}\s{1}\d{1}$",line) is not None:
                                    originalxyz = True
                                elif re.match("^\s*END\n",line) is not None:
                                    originalxyz = False
                                    break
                                if originalxyz:
                                    xyz_original.append(line.strip())
                            xyz_original = list(xyz_original[1:-1])
                            f.seek(0)
                            for line in ff:
                                for xyz in xyz_original:
                                    if xyz in line:
                                        line = "Remove This"
                                if line != "Remove This":
                                    f.write(line)
                            f.truncate()
                            f.seek(0)
                            ff = f.readlines() #needed for replacing later on
                            f.seek(0)
                            for line in ff:
                                if re.match("^\s*\d{1}\s{1}\d{1}$",line) is None:
                                    f.write(line)
                                else:
                                    f.write(line)
                                    for xyz in xyz_reformatted:
                                        f.write(xyz+"\n")
                            f.truncate()
                            subprocess.call(["rm", po_file])
                            subprocess.call(["rm", log_file])
                            subprocess.call(["rm", chk_file])
                            subprocess.call(["qsub", job_file])
                            subprocess.call(["rm", job_outputfile])
                    elif "Memory allocation failure" in job_allocation:
                        subprocess.call(["rm", po_file])
                        subprocess.call(["rm", log_file])
                        subprocess.call(["rm", chk_file])
                        subprocess.call(["qsub", job_file])
                        subprocess.call(["rm", job_outputfile])
                    elif "Maximum iterations reached" in job_allocation:
                        for l in rf:
                            if "Rotational constants (GHZ)" in l or "Cartesian Forces" in l:
                                p = True # results from last iteration begin
                            elif "Number     Number       Type" in l:
                                p = False # results from last iteration end
                                break # don't bother reading further
                            if p:
                                xyz_coords.append(l.strip()) # save results
                        xyz_coords = list(reversed(xyz_coords))
                        xyz_coords = xyz_coords[1:-2] #cutting results
                        for xyz in xyz_coords:
                            xyz = xyz.split()
                            xyz[1] = atom_number_dict.get(xyz[1])
                            del xyz[0]
                            del xyz[1]
                            xyz = "   ".join(xyz)
                            xyz_reformatted.append(xyz)
                        if not os.path.exists("notConv"):
                            os.mkdir("notConv")
                            os.chdir("notConv")
                        else:
                            os.chdir("notConv")
                        with open(f'{molecule_name}_convError.xyz',"w") as f:
                            f.write(f'{len(xyz_coords)}\n\n')
                            for xyz in xyz_reformatted:
                                f.write(f'{xyz}\n')            
                        os.chdir(cwd)
                    elif "Completed Successfully" in job_allocation:
                        for l in rf:
                            if "Rotational constants (GHZ)" in l or "Cartesian Forces" in l:
                                p = True # results from last iteration begin
                            elif "Number     Number              X" in l:
                                p = False # results from last iteration end
                                break # don't bother reading further
                            if p:
                                xyz_coords.append(l.strip()) # save results
                        xyz_coords = list(reversed(xyz_coords))
                        xyz_coords = xyz_coords[1:-2] #cutting results
                        for xyz in xyz_coords:
                            xyz = xyz.split()
                            xyz[1] = atom_number_dict.get(xyz[1])
                            del xyz[0]
                            del xyz[1]
                            xyz = "   ".join(xyz)
                            xyz_reformatted.append(xyz)
                        if not os.path.exists("Conv"):
                            os.mkdir("Conv")
                            os.chdir("Conv")
                        else:
                            os.chdir("Conv")
                        with open(f'{molecule_name}_opt.xyz',"w") as f:
                            f.write(f'{len(xyz_coords)}\n\n')
                            for xyz in xyz_reformatted:
                                f.write(f'{xyz}\n')   
                        os.chdir(cwd) 
            print(job_ID+" "+molecule_output_name+" "+job_allocation)
class OrcaClass:
    def __init__(self):
        cwd = os.getcwd()
        self.cwd = cwd
    def createjob(self):
        with open("job.template", "r") as f:
            job_template = f.read()
        with open("input.template", "r") as f:
            input_template = f.read()
            xyzlist = glob.glob("*xyz")
            for xyz in xyzlist:
                molname = xyz.split(".xyz")[0]
                job_file = job_template.replace("job.template",molname)
                input_file = input_template.replace("job.template",molname)
                if not os.path.exists(molname):       
                    os.mkdir(molname)
                    os.chdir(molname)
                else:
                    os.chdir(molname)                 
                with open(molname + ".job", "w") as f:
                    f.write(job_file)
                with open(molname + ".inp", "w") as f:
                    f.write(input_file)
                os.chdir(self.cwd)
                shutil.copyfile(xyz, f'{molname}/{xyz}')
    def calcenergy(self):
        log_file_list = []
        energies_list = []
        for dirpath, dirnames, filenames in os.walk("."):
            for filename in [f for f in filenames if f.endswith(".out")]:
                log_file_list.append(os.path.join(dirpath, filename))
        if not "orcaenergies.csv" in glob.glob("*.csv"):
            pass
        else:
            os.remove("orcaenergies.csv")
        for log_file in log_file_list:
            totalenergy = "Not Found"
            refenergy = "Not Found"
            ccsdenergy = "Not Found"
            tenergy = "Not Found"
            with open(log_file, "r") as f:
                log_lines = f.readlines()
                for line in log_lines:
                    if "E(CCSD(T))" in line:
                        totalenergy = line.split()[2]
                    if "E(0)" in line:
                        refenergy = line.split()[2]
                    if "E(CORR)(corrected)" in line:
                        ccsdenergy = line.split()[2]
                    if "Triples Correction (T)" in line:
                        tenergy = line.split()[4]
                    correnergy = ccsdenergy + tenergy 
                energies_list.append(f'{log_file},{totalenergy},{refenergy},{ccsdenergy},{tenergy},{correnergy}')
        with open("orcaenergies.csv", "w") as f:
            f.write("Molecule,Total Energy,Reference Energy,CCSD Energy,T Energy,Correlation Energy\n")
            for energy in energies_list:
                f.write(f'{energy}\n')        
class GAMESSClass:
    pass
class MiscFunctions:
    def __init__(self):
        cwd = os.getcwd()
        self.cwd = cwd
    def submitjobs(self):
        print("Hello")
        print("System Used:\nSlurm = 1\nPBS = 2")
        system = int(input())
        print("How many subdirs?")
        subdirs= int(input())
        for root,dirs,files in os.walk(self.cwd):
            if root[len(self.cwd):].count(os.sep) < subdirs:
                for f in files:
                    if f.endswith(".job"):
                        if system == 1:
                            os.chdir(root)
                            subprocess.call(['sbatch','-A','cop@cpu',f])
                            os.chdir(self.cwd)
                        elif system == 2:
                            os.chdir(root)
                            subprocess.call(['qsub',f])
                            os.chdir(self.cwd)
software_map = {"1": GaussianClass,"2": Psi4Class, "3": OrcaClass,"4": CrestClass, "5": GAMESSClass, "6": MiscFunctions}
if program_type in software_map:
    pass
else: 
    print("Software not Found")
    quit()
software = software_map[program_type]()
if program_type == "4":
    software.unpack()
elif program_type == "1":
    print("Function to Perform:\nCreate jobs = 1\nSort Jobs = 2")
    function = input()
    if function == "1":
        software.createjob()
    elif function == "2":
        software.processjob()
elif program_type == "2":
    print("Function to Perform:\nCreate Jobs = 1\nCalculate Energy = 2")
    function = input()
    if function == "1":
        software.createjob()
    elif function == "2":
        software.calcenergy()
elif program_type == "3":
    print("Function to Perform:\nCreate Jobs = 1\nCalculate Energy = 2")
    function = input()
    if function == "1":
        software.createjob()
    elif function == "2":
        software.calcenergy()
elif program_type == "6":
    print("Function to Perform:\nSubmit Jobs = 1")
    function = input()
    if function == "1":
        software.submitjobs()
    if function == "2":
        software.test()
