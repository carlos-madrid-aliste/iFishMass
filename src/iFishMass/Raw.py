from fileinput import filename
import os
import glob
import numpy as np
import logging

log = logging.getLogger(__name__)
    
class Raw:
    def __init__(self, location, debug=False) -> None:
        """ Object initialization.
            
            Parameters:
            ----------
            location: directory path 
            debug:  optional parameter (boolean) for debbuging purposes.
                    set to False for default.
            
            Return: a Raw object.
        """
        self.location = location
        self.subdirs = set()
        self.data = []
        self.debug = debug

        # find sub-directories containing CSV files only and add it to subdirs set.
        # iterate directory
        for path in os.listdir(self.location):
            self.debug and print(f"path={path}, dir={os.path.join(self.location,path)}")
            #log.debug("path={}, dir={}".format(path, os.path.join(self.location,path)))
            
            for (root,dirs,file) in os.walk(os.path.join(self.location,path)):
                for f in file:
                    self.debug and print(f"\troot={root} dirs={dirs} file={f}")
                    # check only for CSV files and discard eveything else.
                    if f.endswith('.csv') :
                        self.subdirs.add(root)
        
        # We need to add some error checking (throw an exception) in case of:
        # 1. A non-existing directory 
        # 2. or directory with no CSV files is used as argument.

    def list_csv_files(self, dir):
        """ List all CSV files in a given directory (dir)
            
            Parameters:
            -----------
            dir: dir_path (string)

            Return:
            list of csv files(list)
        """
        csv_files = list()
        for (root,dirs,file) in os.walk(os.path.join(self.location, dir)):
            self.debug and print(f"\troot={type(root)} dirs={type(str(dirs))} file={type(str(file))}")

            for f in file:
                if f.endswith('.csv') or f.endswith('.CSV'):
                    f = os.path.join(str(root), str(f))
                    csv_files.append(f) 
        return csv_files
        
    def __str__(self):
        """ string representation of the RAW object.
            self.data is not included in the printing. 
        """
        return f"location={self.location}, subdirs={self.subdirs}, debug={self.debug}"
    
    def save_to_csv(self, output_filename, header=True):
        """ Save content of self.data (list of lists) to a text file in CSV format.
            
            Parameters:
            -----------
            output_filename: name of the output filename (string).
                
            header: boolean that determines the writing of the header. Optional argument.
                    header names
                    'M/Z', 'Experimental_M/Z', 'INTENSITY', 'RAW_FILE_NAME', 'IN_FILE (SCAN)'
            
            Throw an exception if ouput file cannot to saved?
        """
        import csv
        
        if len(self.data) == 0:
            assert len(self.data) >= 0, "save_to_csv. data is empty. Do some filtering before saving to csv." 

        #field_names = ['M/Z', 'Experimental_M/Z', 'INTENSITY', 'RAW_FILE_NAME', 'IN_FILE (SCAN)']
        field_names = ['M/Z', 'EXPERIMENTAL_M/Z', 'INTENSITY', 'SAMPLE', 'FILE']
        
        # By default an additional blank line between two subsequent rows is added to the CSV output file.
        # To remove the blank line, you pass the argument newline='' to the open function.
        with open(output_filename, mode='w', newline='') as csv_file:
            csv_file_writer = csv.writer(csv_file, delimiter=',')
            
            header and csv_file_writer.writerow(field_names) # write the header
            csv_file_writer.writerows(self.data)  # write the data
    
    def print_data(self):
        """ print the values of the list self.data 
        """
        import pprint
        pp = pprint.PrettyPrinter(indent=4)
        
        if len(self.data) == 0:
             assert len(self.data) >= 0, "print_data. data is empty. Do some filtering before saving to csv." 
        pp.pprint(self.data)
        
    def get_highest_intensities_per_raw(self, ppm_tolerance, list_of_masses):
        """ get the highest intensity value for a given m/z in all scans of a given RAW file.
        Searching for the highest intensity is done using a tolerance (ppm_tolerance)
        
        A single directory represents a RAW file. Scans associated to a RAW file
        are represented as CSV files inside the directory.
        
        Parameters:
        -----------
        ppm_tolerance  = ppm tolerance value.
        list_of_masses = a list containing a list of m/z values. 
        
        Return:
        list of lists
        [theoretical_mz, experimental_mz, max_intensity, dir, filename])

        [   
            [   6983.36,'595.74756','15172.969','C:\\temp\\MERCK\\CSV\\20210910_Jenny_Merck_Expt1_DI_B4_D6_2','15.csv'],
            [   6983.36,'441.20245','7860.515','C:\\temp\\MERCK\\CSV\\20210910_Jenny_Merck_Expt1_DI_B2_B6_1','170.csv'],
        """
        import numpy as np
        import os
        
        assert ppm_tolerance >= 0, "mz_tolerance must be a positive scalar."
        self.data = []
        
        for mz in list_of_masses:
            self.debug and print(f"looking for mz = {mz}")

            for dir in self.subdirs:
                self.debug and print(f"dir= {dir}") 
                
                my_mzs, my_intensities, my_filenames  = self.filter_by_mz_per_raw(dir, ppm_tolerance, mz)
                
                # error checking
                # skipping directories (RAW) having m/z not within the ppm_tolerance
                if len(my_mzs) == 0: # and len(my_intensities) == 0 and len(my_filenames) == 0:
                    continue

                masses      = np.array(my_mzs).astype(float)  
                intensities = np.array(my_intensities).astype(float)  
                filenames   = np.array(my_filenames).astype(str)
                
                intensities_to_keep = np.array(intensities).astype(float)  
                
                ## get the index of maximum intensity for the massess to keep.
                # get the index of maximum element in numpy array
                result = np.where( intensities_to_keep == np.amax(intensities_to_keep) )
                if self.debug:
                    print('Returned tuple of arrays :', result)
                    print('List of Indices of maximum element :', result[0])   
                #index_max_element = result[0] # it is a tuple
                index_max_element = result[0][0]
                self.debug and print(f"index_max_element={index_max_element}  type={type(index_max_element)}")
                
                # get the filename and experimental m/z
                mz_experimental    = masses[index_max_element]
                max_intensity      = intensities[index_max_element]
                full_path_filename = filenames[index_max_element]
                filename           = os.path.basename(full_path_filename)
            
                self.debug and print(f"SEARCH FOR mz={mz}, mzexp={mz_experimental}, max_inten={max_intensity}, dir={dir}, filename={filename}")
                self.data.append([mz, mz_experimental, max_intensity, dir, filename])
                #self.data.append([mz, mz_experimental, max_intensity, full_path_filename])

        if self.debug:
            print("THIS IS DATA")        
            print(self.data)        
        
        
    def filter_by_mz_per_raw(self, dir, ppm_tolerance, mz):
        """ Filter m/z, within a specified ppm_tolerance, in all scans of a given RAW file.
        
        A single directory represents a RAW file. Scans associated to a RAW file
        are represented as CSV files inside the directory.
            
            Parameters
            ---------
            dir :  directory name (string) representing a single RAW file.
                Inside of directory there are CSV files corresponding to 
                the different spectrum of the RAW file.
            
            ppm_tolerance: ppm_tolerance (integer)
            mz :  m/z value (float)

            Return
            ------
            A list of lists containing all intensity values for a specific m/z within a 
            ppm_tolerance. An empty list is returned if m/z is out of range.
            [
                [intensity0, intensity1 ... intensityn],
                [filename0, filename1 ... filenamen]
                ...
            ]
            
        """
        import numpy as np
        import os
        import pprint
        pp = pprint.PrettyPrinter(indent=4)
        
        assert ppm_tolerance >= 0, "mz_tolerance must be a positive scalar."
        assert mz >= 0, "mass_in must be a positive scalar"

        # listing all scan files.
        my_list = self.list_csv_files(dir) 
        self.debug and print(f"list_csv_files output = {my_list}")
        
        # loading scan files into a bidimensional list.
        filename_mz_intensity_list =  self.load_all_files_in_memory(my_list)
        self.debug and pp.pprint(filename_mz_intensity_list)
        
        # extract m/z, intensities and filenames from the bidimensional list and keep m/z, intensity
        # within a given ppm_tolerance range.
        filenames   = [ my_list[0] for my_list in filename_mz_intensity_list ]
        masses      = [ my_list[1] for my_list in filename_mz_intensity_list ]
        intensities = [ my_list[2] for my_list in filename_mz_intensity_list ]
        
        self.debug and print(f"type={type(masses)} type={type(intensities)}")
        
        # ppm calculation
        masses      = np.array(masses).astype(float)  
        intensities = np.array(intensities).astype(float)  
        filenames   = np.array(filenames).astype(str)
        # lists were converted into a numpy array to avoid this error
        # TypeError: unsupported operand type(s) for -: 'list' and 'float'
        
        masses_to_keep = (np.abs(masses - mz) / mz ) * 1_000_000
        masses_to_keep = masses_to_keep <= ppm_tolerance
        # At this point masses_to_keep is an array of booleans. 
        # The indexes of all masses within the ppm_tolerance range have a True value. 
        
        # Error checking
        if not np.any(masses_to_keep):
            return [], [] , []
        
        # get the mz, intensities and filenames of the masses to keep.
        m_to_keep = masses[masses_to_keep]
        i_to_keep = intensities[masses_to_keep]
        f_to_keep = filenames[masses_to_keep]
        
        return m_to_keep, i_to_keep, f_to_keep
        
    def load_all_files_in_memory(self, list_of_files):
        """ Load a list of CSV files, containing m/z and intensities pairs,
        into a bidimensional list. 
            
            Parameters:
            ----------
            list_of_files: list of CSV files containing mz, intensity (list)
            
            Return
            ------
            A list of list. Inner list contains three elements:
                first : csv filename
                second: experimental m/z
                third:  experimental intensity

                A CSV file might contains multiple m/z, intensity pairs in that case
                file name is repeated in the inner list first field. See below for details.
                [   
                    [   '/MERCK/CSV/20210910_Jenny_Merck_Expt1_DI_B4_D6_1/102.csv',
                        '1190.5034',
                        '15870.921'],
                    [   '/MERCK/CSV/20210910_Jenny_Merck_Expt1_DI_B4_D6_1/102.csv',
                        '1180.5034',
                        '15970.921'],    
                    [   '/MERCK/CSV/20210910_Jenny_Merck_Expt1_DI_B4_D6_1/2.csv',
                        '441.2022',
                        '7151.1636'],
                    [   '/MERCK/CSV/20210910_Jenny_Merck_Expt1_DI_B4_D6_1/308.csv',
                        '441.2029',
                        '10327.367'],
                    [   '/MERCK/CSV/20210910_Jenny_Merck_Expt1_DI_B4_D6_1/32.csv',
                        '587.93494',
                        '4325.3374'],
                    [   '/MERCK/CSV/20210910_Jenny_Merck_Expt1_DI_B4_D6_1/32.csv',
                        '595.7479',
                        '8305.278'],
                    ....
                ] 
        """
        import csv
        import pprint
        pp = pprint.PrettyPrinter(indent=4)
        
        self.debug and print(f"load_all_files_in_memory")
        two_d_list = list()
        
        for file in list_of_files:
            self.debug and print(f"opening={file}")
            with open(file, 'r') as fh:
                reader = csv.reader(fh)
                next(reader) # skip the header
                
                for row in reader:
                    my_list = list()
                    experimental_mz = row[0]
                    intensity = row[1]

                    my_list.append(file)
                    my_list.append(experimental_mz)
                    my_list.append(intensity)
                     
                    two_d_list.append(my_list)
        return two_d_list
    
    def list_is_empty(self, my_list):
        """ Check if a bidimensional python list it contains empty lists.
        
        Parameters:
        ----------
        bidimensional_list: list of lists
        
        Return
        ------
        return True if list of list has empty array elements, otherwise False.
        
        my_list = [ [], [], [] ] 
        if is_empty(my_list):
            print('Array IS empty')
        """
        for inner_list in my_list:
            self.debug and print(f"inner_list={inner_list}")
            if len(inner_list) >= 1:
                return False
        return True      

    def intensities_among_all_raw_files(self, ppm_tolerance, list_of_masses):
        """ List all intensities for a given list of masses (m/z) in each RAW file.

        Each directory represents a RAW file. Inside of each subdir there is a bunch of CSV files 
        that represent different scans or spectrum for that RAW file.

        ppm_tolerance  = ppm value 
        list_of_masses = list()
        """
        import numpy as np
        import itertools

        self.data = []
        for mz in list_of_masses:
            self.debug and print(f"looking for mz = {mz}")

            for dir in self.subdirs:
                self.debug and print(f"dir= {dir}") 
                my_mzs, my_intensities, my_filenames  = self.filter_by_mz_per_raw(dir, ppm_tolerance, mz)
                
                # iterating over multiples list at a time
                for (my_mz, my_intensity, full_path_filename) in zip(my_mzs, my_intensities, my_filenames):
                    filename = os.path.basename(full_path_filename)
                    #self.data.append([mz, my_mz, my_intensity, full_path_filename, filename])
                    self.data.append([mz, my_mz, my_intensity, dir, filename])
                    
        if self.debug:
            print("THIS IS DATA")        
            print(self.data)        

    def get_highest_intensity_among_all_raw_files(self, ppm_tolerance, list_of_masses):
        """ List the highest intensities for a given list of masses (m/z) among all RAW files.

        Each directory represents a RAW file. 
        Inside of each subdir there is a bunch of CSV files 
        that represent differents scans or spectrum for that RAW file.

        ppm   = ppm
        list_of_masses = list()
        """
        import numpy as np
        import os
        import pprint
        pp = pprint.PrettyPrinter(indent=4)
        
        assert ppm_tolerance >= 0, "mz_tolerance must be a positive scalar."
        self.data = []
        for mz in list_of_masses:
            self.debug and print(f"looking for mz = {mz}")

            all_mzs         = list()    
            all_intensities = list()    
            all_filenames   = list()    
            for dir in self.subdirs:
                self.debug and print(f"dir= {dir}") 
                
                my_mzs, my_intensities, my_filenames  = self.filter_by_mz_per_raw(dir, ppm_tolerance, mz)
                # error checking
                # skipping directories (RAW) having m/z not within the ppm_tolerance
                if len(my_mzs) == 0: #and len(my_intensities) == 0 and len(my_filenames) == 0:
                    continue
                
                for m in my_mzs: 
                    all_mzs.append(m)
                for i in my_intensities: 
                    all_intensities.append(i)
                for f in my_filenames: 
                    all_filenames.append(f)
            
            # skipping directories (RAW) having m/z not within the ppm_tolerance
            if len(all_mzs) == 0: #and len(my_intensities) == 0 and len(my_filenames) == 0:
                continue
                
            self.debug and print("PRINTING ALL_MZS")         and pp.pprint(all_mzs)
            self.debug and print("PRINTING ALL_INTENSIITES") and pp.pprint(all_intensities) 
            self.debug and print("PRINTING ALL_FILENAMES")   and pp.pprint(all_filenames)

            all_mzs         = np.array(all_mzs).astype(float)  
            all_intensities = np.array(all_intensities).astype(float)  
            all_filenames   = np.array(all_filenames).astype(str)
            
            intensities_to_keep = np.array(all_intensities).astype(float)  
            
            ## get the index of maximum intensity for the massess to keep.
            # get the index of maximum element in numpy array
            result = np.where( intensities_to_keep == np.amax(intensities_to_keep) )
            if self.debug:
                print('Returned tuple of arrays :', result)
                print('List of Indices of maximum element :', result[0])   
            #index_max_element = result[0] # it is a tuple
            index_max_element = result[0][0]
            self.debug and print(f"index_max_element={index_max_element}  type={type(index_max_element)}")
            
            # get the filename and experimental m/z
            mz_experimental    = all_mzs[index_max_element]
            max_intensity      = all_intensities[index_max_element]
            full_path_filename = all_filenames[index_max_element]
            filename           = os.path.basename(full_path_filename)
        
            self.debug and print(f"SEARCH FOR mz={mz}, mzexp={mz_experimental}, max_inten={max_intensity}, dir={dir}, filename={filename}")
            #self.data.append([mz, mz_experimental, max_intensity, dir, filename])
            self.data.append([mz, mz_experimental, max_intensity, full_path_filename])

        if self.debug:
            print("THIS IS DATA")        
            print(self.data)

    def reshape_long_to_wide(self, csv_filename):
        """ Take a CSV file in long-format and reshape it to wide-format.
        the list of lists in self.data and build a wide table for printing.
        
        Wide-tables have conditions in column_names and sample_names in the rows. (wide form versus long form)
        """                
        import pandas as pd
        import numpy as np
            
        # read some columns of the csv file
        # original columns 'M/Z', 'Experimental_M/Z', 'INTENSITY', 'SAMPLE'
        hi_per_raw = pd.read_csv(csv_filename, usecols=['M/Z', 'INTENSITY', 'SAMPLE'])
        
        #hi_per_raw.columns = ['M/Z', 'INTENSITY', 'SAMPLE']
        #print("THIS IS PANDAS")
        #print(hi_per_raw.head())

        # reshape from long to wide format
        df_wide = pd.pivot(hi_per_raw, index=['SAMPLE'], columns='M/Z', values='INTENSITY')  
        
        # re-arrange the new columns in the correct order
        cols = hi_per_raw['M/Z'].unique()  
        df_wide = df_wide[cols]

        # filling missing values with zero.
        df_wide.fillna(0, inplace=True)
        
        # write to CSV file
        df_wide.to_csv("wide.csv")     
    
    def long_to_wide(self, csv_filename):
        """ Reshape a list of lists (self.data) to wide-format for printing.
        Wide-tables have conditions in column_names and sample_names in the rows. 
        
        Long-format; each row in the table represents a single observation.
        """                
        import pandas as pd
        import numpy as np
            
        # create the panda dataframe from the list of lists
        df = pd.DataFrame(self.data, columns=['M/Z', 'EXPERIMENTAL_MZ', 'INTENSITY', 'SAMPLE', 'FILE'])
        
        # drop a column from data frame on the original object.
        # inplace=True means the operation would work
        # axis=1 means we are dropping the column, not the row 
        # def df['EXPERIMENTAL_MZ']
        df.drop('EXPERIMENTAL_MZ', inplace=True, axis=1)
        #df.to_csv('all_rows.csv')     

        # reshape from long to wide format
        # ValueError: Index contains duplicate entries, cannot reshape
        #df_wide = pd.pivot_table(df, index='SAMPLE', columns='M/Z', values='INTENSITY', aggfunc=np.mean, fill_value=0)  
        df_wide = pd.pivot(df, index=['SAMPLE'], columns='M/Z', values='INTENSITY')  
        #print(df_wide.head())
        #print(df_wide.columns)
        #print(df_wide.keys())

        # appending -M/Z to numeric column names
        new_names = [ str(col_name) + '-M/Z' for col_name in df_wide.columns]    
        # rename all columns at once
        df_wide.columns = new_names

        # re-arrange the new columns in the correct order
        df2 = df_wide.sort_index(axis=1, ascending=True)

        # filling missing values with zero.
        df2.fillna(0, inplace=True)
        
        # write to CSV file
        df2.to_csv(csv_filename)

if __name__ == '__main__':
    import os
    import time
    import sys
    DEBUG = False
    
    start = time.time()
    start_cpu_time = time.process_time()
    
    #/run/media/local/NEW\ VOLUME/MERCK/CSV/
    #dir_path = os.path.join('C:',os.sep ,'run','media', 'local', 'NEW VOLUME', 'MERCK', 'CSV')
    #dir_path = os.path.join('C:',os.sep ,'temp','MERCK', 'CSV')
    dir_path = os.path.join('C:',os.sep ,'temp','MERCK_10_22_2021')
    #dir_path = os.path.join('C:',os.sep ,'tmp','MERCK_10_22_2021')
    
    print(f"dir_path={dir_path}")
    r1 = Raw(dir_path, debug=DEBUG)
    #r1 = Raw(dir_path, debug=True)
    #print(r1)

    list_of_masses = [881.39739, 1761.78747, 587.93404,441.20236,1189.48879,1190.49607,595.75169,397.50357,1296.68481,648.84607,432.89982,324.92669]
    ppm_tolerance = 10

    #list all intensities for a given list of masses (m/z) in each RAW file.
    r1.intensities_among_all_raw_files(ppm_tolerance, list_of_masses)
    output_filename = "intensities_among_all_raw.csv"
    r1.save_to_csv(output_filename)
    DEBUG and print(f"This is DATA {r1.print_data()}")
    print(f"{output_filename} saved!")

    r1.get_highest_intensity_among_all_raw_files(ppm_tolerance, list_of_masses)
    output_filename = "highest_intensities_among_all_raw.csv"
    r1.save_to_csv(output_filename)
    DEBUG and print(f"This is DATA {r1.print_data()}")
    print(f"{output_filename} saved!")
    
    end = time.time()
    elapsed_time = end - start
    process_time = end - start_cpu_time
    print(f"Running time={elapsed_time} seconds")
    print(f'Wall time: {time.strftime("%H:%M:%S", time.gmtime(elapsed_time))}')
    print(f'CPU time: {time.strftime("%H:%M:%S", time.gmtime(process_time))}')
