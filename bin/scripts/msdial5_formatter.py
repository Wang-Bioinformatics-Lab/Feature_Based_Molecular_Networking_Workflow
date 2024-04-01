#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 3/24/2024

@author: mingxun wang
@purpose: to convert the MS-DIAL5 file into a compatible format
"""
import pandas as pd
import sys
from pyteomics import mgf

def convert_to_feature_csv(input_filename, output_filename):
    input_format = pd.read_csv(input_filename, sep='\t', skiprows=4)

    #Check IMS data columns and drop them
    if 'Average drift time' in input_format.columns:
        input_format = input_format.drop(['Average drift time','Average CCS'], axis=1)

    #Continue with the processing
    headers = list(input_format.keys())
    sample_names = headers[22:]

    input_records = input_format.to_dict(orient="records")
    output_records = []

    for record in input_records:
        scan = int(record["Alignment ID"]) + 1
        mz = record["Average Mz"]
        rt = record["Average Rt(min)"]

        output_record = {}
        output_record["row ID"] = str(scan)
        output_record["row m/z"] = str(mz)
        output_record["row retention time"] = str(rt)

        for sample_name in sample_names:
            output_record[sample_name + " Peak area"] = record[sample_name]

        output_records.append(output_record)

    output_headers = ["row ID", "row m/z", "row retention time"]
    output_headers += [sample_name + " Peak area" for sample_name in sample_names]

    output_df = pd.DataFrame(output_records)
    output_df.to_csv(output_filename, sep=",", index=False, columns=output_headers)

    return

def convert_mgf(input_mgf, output_mgf):
    # Write the output MGF file
    with open(output_mgf, "w") as out:
        
        # Cleaning mgf file
        temp_filename = "temp.mgf"

        with open(temp_filename, "w") as outtemp:
            for line in open(input_mgf):
                if "Num Peaks" in line:
                    continue

                outtemp.write(line)

        # reading from pyteomics
        scan = 0
        with mgf.read(temp_filename) as reader:
            for spectrum in reader:
                
                scan += 1
                out.write("BEGIN IONS\n".format(scan))
                out.write("TITLE=SCAN={}\n".format(scan))
                out.write("SCANS={}\n".format(scan))
                # writing out precursor m/z
                out.write("PEPMASS={}\n".format(spectrum["params"]["pepmass"][0]))

                # Figuring out the max intensity and keeping only things above 1% base peak
                try:
                    max_intensity = max(spectrum["intensity array"])

                    for mz, intensity in zip(spectrum["m/z array"], spectrum["intensity array"]):
                        ratio = intensity / max_intensity
                        if ratio > 0.01:
                            out.write("{:.4f} {:.4f}\n".format(mz, intensity))
                except:
                    pass

                out.write("END IONS\n\n".format(scan))

            # # reading from file
            #


        #for line in open(input_mgf):
            # if "TITLE" in line:
            #     scan += 1
            #     out.write("TITLE=SCAN={}\n".format(scan))
            #     out.write("SCANS={}\n".format(scan))
            #     continue
                
            # if "ION" in line:
            #     continue
            # if "CHARGE" in line:
            #     continue
            # if "RTINMINUTES" in line:
            #     continue
            # if "Num Peaks" in line:
            #     continue

            # out.write(line)



    return
    

if __name__=="__main__":
    # there should be obly one input file
    convert_to_feature_csv(sys.argv[1], sys.argv[2])
