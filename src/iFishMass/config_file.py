class config_file:
    def __init__(self, location, debug=False) -> None:
        """ Object initialization.
            
            Parameters:
            ----------
            location: directory path 
            debug:  optional parameter (boolean) for debbuging purposes.
                    set to False for default.
            
            Return: a config_file object.
        """
        self.location = location
        self.debug = debug
        # return error if INI files doe snot exists
        # return error if INI-like format is not provided

    def read_ini(self):
        """ Read INI file and return a dictionary with the key values.
            Parameters:
            -----------
            ini_file: 
                file path to INI file
            Return:
        """
        import sys
        import configparser
        
        # reading INI file
        config = configparser.ConfigParser()
        ini_file = self.location
        config.read(ini_file)
        
        print(f"Reading INI file={ini_file}")
        config_values = dict() 
        try:
            data_folder = config.get('data_folder', 'location')
            output = config.get('output', 'location')
            ms_level = config.get('ms_level', 'level')
            ppm = int(config.get('ppm', 'value'))
            debug = config.getboolean('debug','debug')
            
            # read list of masses from INI file. masses are loaded into
            # a set to remove reduandancy. 
            my_masses = set()
            for mass in config['list_of_masses']:
                #my_masses.append(float(config['list_of_masses'][mass]))
                my_masses.add(float(config['list_of_masses'][mass]))
            
            # convert list into a set to remove redundancy.
            #internal_standard = { x for x in my_list} # set comprehension
            internal_standard = set()
            for mass in config['internal_standard']:
                internal_standard.add(float(config['internal_standard'][mass]))
            
            modified_peptides = set()
            for mass in config['modified_peptides']:
                modified_peptides.add(float(config['modified_peptides'][mass]))
            
            unmodified_peptides = set()
            for mass in config['unmodified_peptides']:
                unmodified_peptides.add(float(config['unmodified_peptides'][mass]))
            
            config_values['data_folder']  = data_folder
            config_values['output'] = output
            config_values['ms_level'] = ms_level
            config_values['ppm']   = ppm
            config_values['debug'] = debug
            config_values['list_of_masses'] = my_masses
            config_values['internal_standard'] = internal_standard
            config_values['modified_peptides'] = modified_peptides
            config_values['unmodified_peptides'] = unmodified_peptides
        except:
            print('Could not read configuration file')
            sys.exit(1)
        return config_values 
    
    def write_demo_ini(self):
        import sys
        import configparser
        
        config = configparser.ConfigParser()
        ini_file = self.location
        
        config.add_section('data_folder')
        config['data_folder']['location']=r"C:\Users\cmadrid\Downloads\Merck_Expt6_DI_LODLOQ_Aug_2022\MZXML"
        
        config.add_section('ms_level')
        config['ms_level']['level']='1'

        config.add_section('ppm')
        config['ppm']['value']='10'
        
        config.add_section('list_of_masses')
        config['list_of_masses']['value1']='881.39739'
        config['list_of_masses']['value2']='1761.78747'
        config['list_of_masses']['value3']='587.93404'
        config['list_of_masses']['value4']='441.20236'
        config['list_of_masses']['value5']='1189.48879'
        config['list_of_masses']['value6']='1190.49607'
        config['list_of_masses']['value7']='595.75169'
        config['list_of_masses']['value8']='397.50357'
        config['list_of_masses']['value9']='298.37951'
        config['list_of_masses']['value10']='1296.68481'
        config['list_of_masses']['value11']='648.84607'
        config['list_of_masses']['value12']='432.89982'
        config['list_of_masses']['value13']='324.92669'
        
        config.add_section('debug')
        config['debug']['debug']="False"

        config.add_section('output')
        config['output']['location']="C:/temp/MERCK_AUGUST"

        config.add_section("internal_standard")
        config['internal_standard']['value1']='1296.68481'
        config['internal_standard']['value2']='648.84607'
        config['internal_standard']['value3']='432.89982'
        config['internal_standard']['value4']='324.92669'
        
        config.add_section("modified_peptides")
        config['modified_peptides']['value1']='881.39739'
        config['modified_peptides']['value2']='1761.78747'
        config['modified_peptides']['value3']='587.93404'
        config['modified_peptides']['value4']='441.20236'

        config.add_section("unmodified_peptides")
        config['unmodified_peptides']['value1']='1189.48879'
        config['unmodified_peptides']['value2']='1190.49607'
        config['unmodified_peptides']['value3']='595.75169'
        config['unmodified_peptides']['value4']='397.50357'
        config['unmodified_peptides']['value5']='298.37951'

        
        try:
            #with open(ini_file) as config_file:
            with open('peak.ini', 'w') as config_file:
                config.write(config_file)    
        except:
            print('Could not write configuration file')
            sys.exit(1)
        
