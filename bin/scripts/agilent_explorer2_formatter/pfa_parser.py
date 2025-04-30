import xml.etree.ElementTree as ET
import zipfile
import csv
from io import BytesIO, TextIOWrapper
import logging

Debug = False

# Function to select, read and validate a .pfa file:
def ReadPFAFile(file_path):

    logging.info("Validating PFA file ...")
    try:
        with zipfile.ZipFile(file_path, 'r') as zip_ref:
            # Print the number of files in the ZipFile
            if Debug:
                logging.info(f"Number of files in the .pfa: {len(zip_ref.namelist())}")

            has_sample_group = False
            has_cef_file = False
            has_MSMS_file = False

            for file_name in zip_ref.namelist():
                if Debug:
                    logging.info(f"{file_name}")
                if file_name.endswith('SampleGroups.tsv'):
                    has_sample_group = True
                if file_name.endswith('.cef'):
                    has_cef_file = True
                if file_name.endswith('CompositeMSMS.msp'):
                    has_MSMS_file = True

    except FileNotFoundError:
        logging.info(f"ERROR: File not found: {file_path}")
        exit()
    except zipfile.BadZipFile:
        logging.info(f"ERROR: Bad zip file: {file_path}")
        exit()
    except Exception as e:
        logging.info(f"ERROR: An error occurred: {e}")
        # exit()
        raise e

    # Validate file_name has SampleGroup.tsv and at least one file ending with '.cef'
    if not has_cef_file:
        logging.info(f"ERROR: {file_path} does not have any compound results (no .cef file(s) found)")
        exit()
    if not has_sample_group:
        logging.info(f"ERROR: {file_path} does not have Sample Grouping file (SampleGroups.tsv file not found)")
        exit()
    if not has_MSMS_file:
        logging.info(f"ERROR: {file_path} does not have a MS/MS spectra (CompositeMSMS.msp file not found) please make sure MS/MS spectra has been extracted can re-export the .pfa.")
        exit()
    if has_sample_group and has_cef_file:
        logging.info(f"        ... {file_path} is valid.")
    return file_path

# Function to read sample grouping information from the .pfa file into memory:
def GetSampleGroupingFromPFA(file_path, Debug=False):
    logging.info("Reading Sample Groups ...")

    SampleGroups = None
    with zipfile.ZipFile(file_path, 'r') as zip_ref:
        for file_name in zip_ref.namelist():
            # open each file from the .pfa file
            with zip_ref.open(file_name) as file:
                # read sample grouping file (structured as a tab separated text file) into memory:
                if file_name == "SampleGroups.tsv":
                    SampleGroups_content = TextIOWrapper(BytesIO(file.read()), encoding='utf-8')
                    SampleGroups_reader = csv.DictReader(SampleGroups_content, delimiter='\t')
                    SampleGroups = list(SampleGroups_reader)

    if Debug and SampleGroups:
        for row in SampleGroups:
            logging.info(row)

    # Check if SampleGroups is not empty
    if SampleGroups:
        logging.info("        ... Sample Groups extracted")
    else:
        logging.info("        ... No Sample Group information found!")
        return None
    
    return SampleGroups

#Function to read the compound results from each sample result file (.cef) in the .pfa file into memory:
def GetCompoundResultsFromPFA(file_path, output_dir):
    logging.info("Extracting Compound Group results ...")
    
    # Initialize object for .pfa contents
    file_name_to_path = {}

    # open .pfa file (a .zip structured file)
    with zipfile.ZipFile(file_path, 'r') as zip_ref:
        total_files = len(list(zip_ref.namelist()))
        for i, file_name in enumerate(zip_ref.namelist()):
            
            # open each file from the .pfa file
            with zip_ref.open(file_name) as file:

                # read each .cef file (compound exchange format, structured as a .xml file) into memory:
                if file_name.endswith('.cef'):
                    if file_name in file_name_to_path:
                        raise ValueError(f"Duplicate filename detected: {file_name}")
                    
                    output_path = f"{output_dir}/{file_name}"
                    with open(output_path, 'wb') as output_file:
                        output_file.write(file.read())
                    
                    file_name_to_path[file_name] = output_path
        
            logging.info(f"        ... Extracted file {i+1} of {total_files}: {file_name}")
    
    # Check if any .cef files were extracted
    if file_name_to_path:
        logging.info("        ... Compound Groups extracted and written to files")
    else:
        logging.info("        ... No Compound Group results found!")
        exit()
    
    return file_name_to_path

# Function to read Composite MS/MS information from the .pfa file into memory:
def GetSCompositeMSMSFromPFA(file_path):
    logging.info("Reading Composite MS/MS ...")

    mspSpectra = None
    with zipfile.ZipFile(file_path, 'r') as zip_ref:
        for file_name in zip_ref.namelist():
            # open each file from the .pfa file
            if file_name == "CompositeMSMS.msp":
                mspSpectra = []
                with zip_ref.open(file_name) as file:
                    lines = TextIOWrapper(file, encoding='utf-8').readlines()
                    spectrum = {}
                    peaks = []
                    prefix_map = {
                        "SYNONYM:$:06": "instrument",
                        "SYNONYM:$:00": "type",
                        "SYNONYM:$:04": "precursormz",
                        "SYNONYM:$:10": "ionsource",
                        "Ion_mode:": "ion_mode",
                        "Charge:": "charge",
                        "Precursor_type:": "precursor_type",
                        "Num Peaks:": "num_peaks"
                    }
                    for line in lines:
                        line = line.strip()
                        if line.startswith("name:"):
                            if spectrum:
                                spectrum['peaks'] = peaks
                                mspSpectra.append(spectrum)
                            name_details = line.split(": ")[1]
                            name_parts = name_details.split(", ")
                            name_value = name_parts[0]
                            if name_value.startswith("Cpd Grp"):
                                name_value = int(name_value.split(" ")[2])
                            spectrum = {
                                'num': name_value,
                                'M': float(name_parts[1].split(" = ")[1]),
                                'RT': float(name_parts[2].split(" = ")[1]),
                            }
                            peaks = []
                        else:
                            for prefix, key in prefix_map.items():
                                if line.startswith(prefix):
                                    value = line.split(": ")[1] if ": " in line else line.split(":")[2][2:]
                                    if key in ["charge", "num_peaks"]:
                                        spectrum[key] = int(value)
                                    elif key == "precursormz":
                                        spectrum[key] = float(value)
                                    else:
                                        spectrum[key] = value
                                    break
                            else:
                                if line:
                                    mz, intensity = map(float, line.replace(';', '').split())
                                    peaks.append((mz, intensity))
                    if spectrum:
                        spectrum['peaks'] = peaks
                        mspSpectra.append(spectrum)
    if mspSpectra:
        logging.info("        ... MS/MS spectra extracted")
    else:
        logging.info("        ... No MS/MS spectra found!")
        exit()
    return mspSpectra

# Function to gather MS1 ion information to be exported, then store information as a list in a tab-seperated text file:
def ExportMS1Ions(PFACompoundResults, output_file):
    
    # Initialize an empty list to store Compound data
    DataToExport = []
    logging.info("Exporting MS1 ions...")
    
    # Iterate through the Compound elements in each tree (there is a XML tree for each file_name)
    for file_name, tree in PFACompoundResults.items():

        # Find the required information needed to in the Export
        for Compound in tree.findall('.//Compound'):
            # Compound information
            CompoundResult = {"File_Name": file_name}
            CompoundResult["num"] = int(Compound.get('num'))
            CompoundResult["Mass"] = float(Compound.find('.//Location').get('m'))
            CompoundResult["RT"] = float(Compound.find('.//Location').get('rt'))
            
            # MS1 Spectra information
            for Spectrum in Compound.findall('.//Spectrum'):
                type = Spectrum.get('type')
                if type == 'FbF' or type == 'MFE':
                    CompoundResult["polarity"] = Spectrum.find('.//MSDetails').get('p')
                    CompoundResult["mz"] = []
                    CompoundResult["intensity"] = []
                    CompoundResult["z"] = []
                    CompoundResult["species"] = []

                    # FbF or MFE MS1 Spectra
                    for MSPeak in Spectrum.findall('.//MSPeaks/p'):
                        CompoundResult["mz"].append(float(MSPeak.attrib['x']))
                        CompoundResult["intensity"].append(float(MSPeak.attrib['y']))
                        CompoundResult["z"].append(float(MSPeak.attrib['z']))
                        CompoundResult["species"].append(MSPeak.attrib['s'])
                    DataToExport.append(CompoundResult)
    if DataToExport:
        logging.info("        ... MS1 ions gathered")
    else:
        logging.info("        ... No FbF or MFE MS1 spectra found!")
        exit()

    # Write the data to a tab-separated text file
    with open(output_file, 'w', newline='') as file:
        writer = csv.writer(file, delimiter='\t')
        
        # Write the header
        writer.writerow(["File_Name", 
                         "Compound Group", 
                         "Mass", 
                         "RT", 
                         "polarity", 
                         "z", 
                         "species",
                         "mz", 
                         "intensity"
        ])
        
        # Write the data
        for row in DataToExport:
            for i in range(len(row["mz"])):
                writer.writerow([
                    row["File_Name"],
                    row["num"],
                    row["Mass"],
                    row["RT"],
                    row["polarity"],
                    row["z"][i],
                    row["species"][i],
                    row["mz"][i],
                    row["intensity"][i]
                ])
    logging.info(f"        ... MS1 Ions saved to {output_file}")

def ExportCompoundResults(PFACompoundResults, ResponseType='Height'):
    
    # Initialize an empty list to store Compound data
    DataToExport = []
    logging.info("Exporting integration results...")
    
    # Iterate through the Compound elements in each tree (there is a XML tree for each file_name)
    for file_name, file_path in PFACompoundResults.items():
        with open(file_path, 'r') as file:
            Sample = file.read()
            tree = ET.ElementTree(ET.fromstring(Sample)).getroot()

            # Find the required information needed to in the Export
            for Compound in tree.findall('.//Compound'):
                # Compound information
                CompoundResult = {"File_Name": file_name}
                CompoundResult["num"] = int(Compound.get('num'))
                CompoundResult["Mass"] = float(Compound.find('.//Location').get('m'))
                CompoundResult["RT"] = float(Compound.find('.//Location').get('rt'))
                
                if ResponseType == 'Height':
                    CompoundResult["Response"] = float(Compound.find('.//Location').get('y'))
                elif Compound.get('algo') == 'FindByMolecularFeature':
                    CompoundResult["Response"] = float(Compound.find('.//Location').get('v'))
                elif Compound.get('algo') == 'FindByFormula':
                    CompoundResult["Response"] = float(Compound.find('.//Location').get('a'))

                DataToExport.append(CompoundResult)
    if DataToExport:
        logging.info("        ... Integration results gathered")
    else:
        logging.info("        ... No integration results found!")
        exit()

    return DataToExport

# Function to gather MS1 ion information to be exported, then store information as a list in a tab-seperated text file:
def ExportCompoundResultsToFile(PFACompoundResults, output_file, ResponseType='Height'):
    
    # Initialize an empty list to store Compound data
    DataToExport = []
    logging.info("Exporting integration results...")
    
    # Iterate through the Compound elements in each tree (there is a XML tree for each file_name)
    for file_name, tree in PFACompoundResults.items():

        # Find the required information needed to in the Export
        for Compound in tree.findall('.//Compound'):
            # Compound information
            CompoundResult = {"File_Name": file_name}
            CompoundResult["num"] = int(Compound.get('num'))
            CompoundResult["Mass"] = float(Compound.find('.//Location').get('m'))
            CompoundResult["RT"] = float(Compound.find('.//Location').get('rt'))
            
            if ResponseType == 'Height':
                CompoundResult["Response"] = float(Compound.find('.//Location').get('y'))
            elif Compound.get('algo') == 'FindByMolecularFeature':
                CompoundResult["Response"] = float(Compound.find('.//Location').get('v'))
            elif Compound.get('algo') == 'FindByFormula':
                CompoundResult["Response"] = float(Compound.find('.//Location').get('a'))

            DataToExport.append(CompoundResult)
    if DataToExport:
        logging.info("        ... Integration results gathered")
    else:
        logging.info("        ... No integration results found!")
        exit()

    # Write the data to a tab-separated text file
    with open(output_file, 'w', newline='') as file:
        writer = csv.writer(file, delimiter='\t')
        
        # Write the header
        writer.writerow(["File_Name", 
                         "Compound Group", 
                         "Mass", 
                         "RT", 
                         "Response"
        ])
        
        # Write the data
        for row in DataToExport:
            writer.writerow([
                row["File_Name"],
                row["num"],
                row["Mass"],
                row["RT"],
                row["Response"]
                ])
    logging.info(f"        ... Integration results saved to {output_file}")