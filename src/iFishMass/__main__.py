""" iFishMass.py 
Filter all XML files by list_of_masses with a specific ppm_tolerance
save the resulting filtered files in CSV format. One file per scan.

python .\peakEXtractor.py --inifile C:/Users/cmadrid/Downloads/peakEXtractor/peak.ini
Author: Carlos Madrid-Aliste 
Date  : August 23, 2022
email : carlos.madrid-aliste@einsteinmed.edu
        creggae@gmail.com
"""
import sys
#import os
#from xml.etree.ElementTree import PI
#from pyteomics import mzxml, auxiliary
#from pyteomics.auxiliary import cvquery
#import spectrum_utils.spectrum as sus
#from Raw import Raw 

def keep_peaks_around_mass_cma(spectrum_in, mass_in, ppm_tolerance):
    import numpy as np
    import copy
    """ Keep peaks that are within mz_tolerance (in ppm) of the mass_in.
        
        Paramaters:
        ----------
        spectrum_in: 
            Input spectrum
        mass_in: 
            mass to filter the peaks
        ppm_tolerance: 
            tolerance of mz values (in ppm)
    """
    if spectrum_in is None:
        return None

    assert ppm_tolerance >= 0, "mz_tolerance must be a positive scalar." 
    assert mass_in >= 0, "mass_in must be a positive scalar"

    spectrum = copy.deepcopy(spectrum_in)   

    mzs, intensities = spectrum['m/z array'], spectrum['intensity array']
    
    # ppm calculation
    peaks_to_keep = (np.abs(mzs - mass_in) / mass_in ) * 1_000_000
    
    peaks_to_keep = peaks_to_keep <= ppm_tolerance
    
    #condition = np.logical_and(spectrum.peaks.mz <= spectrum.peaks.mz,)
    new_mzs, new_intensities = mzs[peaks_to_keep], intensities[peaks_to_keep]
    spectrum['m/z array'] = new_mzs
    spectrum['intensity array'] = new_intensities
    return spectrum

def filter_peaks(spectrum_in, list_of_masses, ppm_tolerance, debug=False):
    import numpy as np
    import copy
    """ Keep peaks that are within mz_tolerance (in ppm) of the list_of_masses
        
        Paramaters:
        ----------
        spectrum_in: 
            Input spectrum
        list_of_masses: 
            list of masses to filter the peaks
        ppm_tolerance: 
            tolerance of mz values (in ppm)
        Return:
            spectrum object (dictionary)    
    """
    if spectrum_in is None:
        return None

    assert ppm_tolerance >= 0, "mz_tolerance must be a positive scalar." 
    assert len(list_of_masses) >= 0, "list_of_masses is empty"

    spectrum = copy.deepcopy(spectrum_in)   
    #mzs, intensities = np.array(spectrum['m/z array']), np.array(spectrum['intensity array'])
    mzs, intensities = spectrum['m/z array'], spectrum['intensity array']
    #if np.all(mzs):
    #    print(f"mzs all elements are True")
    #else:
    #    print(f"mzs all elements are False")

    # list where to store the peaks, initialized with zeros.
    peaks_to_keep = np.zeros(len(mzs), dtype=bool)
    #peaks_to_keep = np.full(len(mzs), False)
    
    for mass_in in list_of_masses:
        debug and print (f"LOOKING for {mass_in} type={type(mass_in)}")
        # ppm calculation
        tmp_peaks_to_keep = (np.abs(mzs - mass_in) / mass_in) * 1_000_000
        tmp_peaks_to_keep = tmp_peaks_to_keep <= ppm_tolerance

        peaks_to_keep = np.logical_or(peaks_to_keep, tmp_peaks_to_keep)
        if debug and np.any(tmp_peaks_to_keep):
            print(f"TRUE values")

    #print(f"PEAKS_TO_KEEP={peaks_to_keep}")
    new_mzs, new_intensities = mzs[peaks_to_keep], intensities[peaks_to_keep]
    spectrum['m/z array'] = new_mzs
    spectrum['intensity array'] = new_intensities
    return spectrum

def spectrum_is_empty(spectrum_in):
    """ Check if spectrum_in is empty
        
        Paramaters:
        ----------
        spectrum_in: 
            Input spectrum
    """
    if spectrum_in is None:
        return 1
    else:
        # spectrum_in is a dictionary. Check that 'm/z array' and 'intensity array'
        # are not empty.
        if len(spectrum_in['m/z array']) == 0 or len(spectrum_in['intensity array']) == 0:
            return 1
        else:
            return 0       
        
def merge(list1, list2):
    """" merge into a list a list of tuples
    using zip() method to merge two list elements and typecasting
    into tuple.
        Parameters:
        ----------
        list1: first list
        list2: second list
        
        Return:
        list
    """
    merged_list = tuple(zip(list1, list2))
    return merged_list

def save_as_csv(spectrum_in, filename):
    import csv
    """ Keep peaks that are within mz_tolerance (in ppm) of the mass_in.
        
        Paramaters:
        ----------
        spectrum_in: 
            Input spectrum
        mass_in: 
            mass to filter the peaks
        ppm_tolerance: 
            tolerance of mz values (in ppm)
    """
    if spectrum_in is None:
        return None

    mzs, intensities = spectrum_in['m/z array'], spectrum_in['intensity array']
    # create a list of tuples. It is easier to write to CSV in this way.
    data = merge(mzs, intensities)
    header = ['m/z', 'intensity']

    with open(filename, 'w', encoding='UTF8', newline='') as f:
        writer  = csv.writer(f)
        writer.writerow(header) # write the header
        writer.writerows(data)  # write multiple rows
    return 1    

def get_mzxml_files(dir_path):
    # list to store mxXML files
    my_list = []
    
    # iterate directory
    for file in os.listdir(dir_path):
        # check if current path is a file
        if os.path.isfile(os.path.join(dir_path, file)):
            my_file = os.path.join(dir_path, file)
            if my_file.endswith('.mzXML'):
                my_list.append(my_file)
    return my_list 

def get_mzxml_files_yield(dir_path):
    import os
    # using a generator expression will buid a generator function.
    # this approach is memory efficient.

    # iterate directory
    for file in os.listdir(dir_path):
        # check only mzXML files
        if os.path.isfile(os.path.join(dir_path, file)):
            my_file = os.path.join(dir_path, file)
            if my_file.endswith('.mzXML'):
                yield my_file

def filter_files(*, input_dir, output_dir, ms_level, ppm_tolerance, debug, list_of_masses):
    """ Filter all XML files by list_of_masses with a specific ppm_tolerance
        save the resulting filtered files in CSV format. One file per scan.

        Parameters:
        -----------
        conf_file: 
            file path to INI file
        Return:
    """
    import sys
    import os
    from tqdm import tqdm
    from pyteomics import mgf, auxiliary
    from xml.etree.ElementTree import PI
    from pyteomics import mzxml, auxiliary
    from pyteomics.auxiliary import cvquery
    import spectrum_utils.spectrum as sus
    import pprint
    pp = pprint.PrettyPrinter(indent=4)

    #for file_xml in get_mzxml_files_yield(input_dir):
    # Wrapping tqdm around an iterable.
    pbar = tqdm(get_mzxml_files_yield(input_dir))
    for file_xml in pbar:
        debug and print(f"XML_FILE={file_xml}")
        
        filename = os.path.basename(file_xml)
        #pbar.set_description("Processing %s" % filename)
        pbar.set_description("Processing %s" % file_xml)
        
        # iterate over all scans, filter and write them
        # into a CSV formmated file.
        with mzxml.read(file_xml) as reader:
            #debug and auxiliary.print_tree(next(reader))
            for spectrum in reader:
                spectrum_ms_level = str(spectrum['msLevel'])
                
                if debug:
                    print(f"LOOP spectrum type={type(spectrum)}")
                    print(f"keys={spectrum.keys()}")
                    debug and print(f"THIS IS spectrum")
                    debug and pp.pprint(spectrum)
                    
                    print(f"msLevel  ={spectrum['msLevel']}")
                    print(f"polarity ={spectrum['polarity']}")
                    print(f"filterLine= {spectrum['filterLine']}")
                    print(f"num ={spectrum['num']}")
                    print(f"id  ={spectrum['id']}")
                    print(f"m/z array ={spectrum['m/z array']}")
                    print(f"intensity array ={spectrum['intensity array']}")
                    print(f"spectrum_ms_level={type(spectrum_ms_level)}")
                    print(f"ms_level_config={type(ms_level)}")
                
                # look for peaks only in the ms_level set in the INI file
                if spectrum_ms_level != ms_level:
                    debug and print(f"SKIPPING spectrum_ms_level={spectrum_ms_level}  ms_level_config={ms_level}")
                    continue

                debug and print(f"list_of_masses={list_of_masses}")
                sp = filter_peaks(spectrum, list_of_masses, ppm_tolerance, debug=debug)
                
                debug and print (f"This is sp = {sp}")
                if spectrum_is_empty(sp):
                    debug and print("SPECTRUM IS EMPTY")
                    continue
                
                debug and print(f"SPECTRUM is not empty {spectrum['num']} {sp['num']}")

                # store csv inside of CSV directory

                # 20210910_Jenny_Merck_Expt1_DI_B4_D6_2.mzXML"
                # remove .mzXML from every mzXML file a create a directory.
                # That directory (20210910_Jenny_Merck_Expt1_DI_B4_D6_2) it represents a single raw 
                # file.
                file_tmp_path = output_dir
                
                # os.path.join(os.getcwd(), 'new_folder', 'file.txt')
                if not os.path.exists(file_tmp_path):
                    os.mkdir(file_tmp_path)
                    print(f"Directory {file_tmp_path} created")
                    
                # removing filename extension (.mzXML) to create a dir name based on filename
                filename = os.path.basename(file_xml)
                filename_without_extension = os.path.splitext(filename)[0]
                debug and print(f"filename_without_extension ={filename_without_extension}")

                # create directory where to store the CSV files
                csv_dir_name = os.path.join(file_tmp_path, filename_without_extension)
                if not os.path.exists(csv_dir_name):
                    os.mkdir(csv_dir_name)
                    debug and print(f"Directory {csv_dir_name} created")

                filename = f"{sp['num']}.csv"
                csv_full_path_name = os.path.join(csv_dir_name, filename)
                debug and print(f"filename = {filename}, full_path={csv_full_path_name}")
                #save_as_csv(spectrum, filename) 
                save_as_csv(sp, csv_full_path_name) 

def read_options(args=sys.argv[1:]):
    import argparse

    # construct the argument parser
    parser = argparse.ArgumentParser(
        description="Easy way to extract M/Z and intensities from mzXML files.",
        epilog='Written by Carlos Madrid-Aliste.'
    )
    
    group = parser.add_mutually_exclusive_group()

    # add the positional arguments to the Argument Parser
    group.add_argument("--inifile", default='peak.ini', help="INI-file location")
    group.add_argument("--printini", action='store_true', help="print a demo peak.ini file to current directory.")
    
    # parse arguments from terminal
    opts = parser.parse_args(args)
    if not opts.inifile:   # manually catching mistakes
        parser.error("INI file is required")
    
    return opts

def dump(*, input_dir, output_dir, ms_level):
    """ Save m/z and intensitites for all scans in a given raw file. Files are
    stored in CSV format and no filtering is performed at all.

        Parameters:
        -----------
        input_dir: dir containing mzxml files.
        output_dir: dir where to store the CSV files.
        ms_level: ms_level (mz or mz/mz)
        Return:
    """
    import sys
    import os
    from tqdm import tqdm
    from pyteomics import mgf, auxiliary
    from xml.etree.ElementTree import PI
    from pyteomics import mzxml, auxiliary
    from pyteomics.auxiliary import cvquery
    import spectrum_utils.spectrum as sus

    #for file_xml in get_mzxml_files_yield(input_dir):
    # Wrapping tqdm around an iterable.
    pbar = tqdm(get_mzxml_files_yield(input_dir))
    for file_xml in pbar:
        debug and print(f"XML_FILE={file_xml}")
        
        filename = os.path.basename(file_xml)
        #pbar.set_description("Processing %s" % filename)
        pbar.set_description("Processing %s" % file_xml)
        
        # iterate over all scans, and write them
        # into a CSV formmated file.
        with mzxml.read(file_xml) as reader:
            #debug and auxiliary.print_tree(next(reader))
            for spectrum in reader:
                spectrum_ms_level = str(spectrum['msLevel'])
                if debug:
                    print(f"LOOP spectrum type={type(spectrum)}")
                    print(spectrum)
                    print(f"msLevel  ={spectrum['msLevel']}")
                    print(f"polarity ={spectrum['polarity']}")
                    print(f"filterLine= {spectrum['filterLine']}")
                    print(f"num ={spectrum['num']}")
                    print(f"id  ={spectrum['id']}")
                    print(f"m/z array ={spectrum['m/z array']}")
                    print(f"intensity array ={spectrum['intensity array']}")
                    print(f"spectrum_ms_level={type(spectrum_ms_level)}")
                    print(f"ms_level_config={type(ms_level)}")
                # look for peaks only in the ms_level set in the INI file
                if spectrum_ms_level != ms_level:
                    debug and print(f"SKIPPING spectrum_ms_level={spectrum_ms_level}  ms_level_config={ms_level}")
                    continue

                # store csv inside of CSV directory

                # 20210910_Jenny_Merck_Expt1_DI_B4_D6_2.mzXML"
                # remove .mzXML from every mzXML file a create a directory.
                # That directory (20210910_Jenny_Merck_Expt1_DI_B4_D6_2) it represents a single raw 
                # file.
                file_tmp_path = output_dir
                
                # os.path.join(os.getcwd(), 'new_folder', 'file.txt')
                if not os.path.exists(file_tmp_path):
                    os.mkdir(file_tmp_path)
                    print(f"Directory {file_tmp_path} created")
                    
                # removing filename extension (.mzXML) to create a dir name based on filename
                filename = os.path.basename(file_xml)
                filename_without_extension = os.path.splitext(filename)[0]
                debug and print(f"filename_without_extension ={filename_without_extension}")

                # create directory where to store the CSV files
                csv_dir_name = os.path.join(file_tmp_path, filename_without_extension)
                if not os.path.exists(csv_dir_name):
                    os.mkdir(csv_dir_name)
                    debug and print(f"Directory {csv_dir_name} created")

                filename = f"{spectrum['num']}.csv"
                csv_full_path_name = os.path.join(csv_dir_name, filename)
                debug and print(f"filename = {filename}, full_path={csv_full_path_name}")
                
                save_as_csv(spectrum, csv_full_path_name) 
                debug and print(f"filename = {filename} SAVED SUCESSFULLY")

"""
def get_mzxml_files_yield(dir_path):
    import os
    # using a generator expression will buid a generator function.
    # this approach is memory efficient.

    # iterate directory
    for file in os.listdir(dir_path):
        # check only mzXML files
        if os.path.isfile(os.path.join(dir_path, file)):
            my_file = os.path.join(dir_path, file)
            if my_file.endswith('.mzXML'):
                yield my_file
"""
def remove_dir_content(my_dir):
    """ remove directory content.
        function used to get rid of the temporary csv files.
    
    Paremeters:
    -----------
        dir: directory containing CSV files.
    """
    import shutil
    import os
    import sys
    print (f"\n=====================================================================")
    print(f"To run properly iFishMass need to empty {my_dir}.")
    print(f"iFishMass will remove the content of {my_dir} now!")
    result = input("Should I proceed [y/n]? ")
    if result == 'y' or result == 'Y':
        print(f"Empty {my_dir} content before running ...")
        for filename in os.listdir(my_dir):
            file_path = os.path.join(my_dir, filename)
            try:
                if os.path.isfile(file_path) or os.path.islink(file_path):
                    print(f"Removing file {file_path} ...")
                    os.unlink(file_path)
                elif os.path.isdir(file_path):
                    print(f"Removing dir {file_path} ...")
                    shutil.rmtree(file_path)
            except Exception as e:
                print(f"Failed to delete {file_path}. Reason {e}")            
    else:
        print("Not doing anything. Bye!")
        exit(1)    

def main():

    import os.path
    import time
    import sys
    import pprint
    import pkg_resources
    
    import iFishMass
    from iFishMass import Raw as r
    from iFishMass import config_file  as cfg
    from iFishMass import DataAnalysis as da
    
    pp = pprint.PrettyPrinter(indent=4)
    DO_PLOTS = False    
    
    start = time.time()

    # call the function to read the argument values
    opts = read_options(sys.argv[1:])
    
    # write demo_ini file if asked by end-user
    if opts.printini: 
        cfg_obj = cfg.config_file(location='peak.ini')
        cfg_obj.write_demo_ini()
        print(f"peak.ini file written to current directory")
        print(f"Please modify peak.ini according to your needs.")
        print(f"Then run iFishMass")
        sys.exit() 
    
    ini_file = opts.inifile
    if not os.path.exists(ini_file):
        sys.exit(f"{ini_file} does not exists in the current directory\nPlease type iFishMass --help")

    # read the configuration INI file and extract its values.
    cfg_obj = cfg.config_file(location=ini_file, debug=False)
    values = cfg_obj.read_ini()
    if values['debug'] == True:
        pp.pprint(values)

    idir  = values['data_folder']
    odir  = values['output']
    level = values['ms_level']
    ppm   = values['ppm']
    debug = values['debug']
    masses= values['list_of_masses']
    internal_standard   = values['internal_standard']
    modified_peptides   = values['modified_peptides']
    unmodified_peptides = values['unmodified_peptides']
    
    if internal_standard and modified_peptides and unmodified_peptides:
        DO_PLOTS=True
    
    # remove temporary CSV files before running analysis
    # CSV files from previous will distort results.
    remove_dir_content(odir)
    
    filter_files(input_dir=idir, output_dir=odir, 
        ms_level=level, ppm_tolerance=ppm, debug=debug, list_of_masses=masses
    )

    print(f'Generating  CSV reports ...')
    output = odir
    r1 = r.Raw(output)
    
    r1.intensities_among_all_raw_files(ppm_tolerance=ppm, list_of_masses=masses)
    output_filename = "intensities_among_all_raw.csv"
    r1.save_to_csv(output_filename)
    print(f"\treport {output_filename} saved!")
    #r1.print_data()
    
    r1.get_highest_intensities_per_raw(ppm_tolerance=ppm, list_of_masses=masses)
    output_filename = "highest_intensities_per_raw.csv"
    r1.save_to_csv(output_filename)
    print(f"\treport {output_filename} saved!")
    r1.long_to_wide(csv_filename='highest_intensities_per_raw_wide.csv')
    #r1.print_data()

    r1.get_highest_intensity_among_all_raw_files(ppm_tolerance=ppm, list_of_masses=masses)
    output_filename = "highest_intensities_among_all_raw.csv"
    r1.save_to_csv(output_filename)
    print(f"\treport {output_filename} saved!")

    if DO_PLOTS:
        print(f"Generating plots ...")
        
        csv_filename = 'highest_intensities_per_raw_wide.csv'
        output_plot_file = 'analysis_plot.xlsx'
        try:
            tmpl_file = pkg_resources.resource_filename(__name__, 'data/template.xlsx')
            d = da.DataAnalysis(
                location = csv_filename, 
                modified_peptides  = modified_peptides, 
                unmodified_peptides= unmodified_peptides, 
                internal_standard  = internal_standard,
                debug = debug,
                template_file = tmpl_file,
                output = output_plot_file
            )

            d.do_analysis() 
            if os.path.exists(output_plot_file):
                print(f"\t{output_plot_file} done!")  
        
        except AssertionError as error:
            print(error)

    end = time.time()
    elapsed_time = end - start
    print()
    print(f'Wall time: {time.strftime("%H:%M:%S", time.gmtime(elapsed_time))}')

if __name__ == '__main__':
    main()
