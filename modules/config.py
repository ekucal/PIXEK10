############################################################
############################################################
######  Here, are initial parameters for PIXE code     #####
######  data for UO2 background removal                #####
############################################################
############################################################

import json
class InitialParameters():
    def __init__(self, _init_file_loc_, _init_name_, _init_material_, _init_crystal_orientation_, _init_axial_channel_, _init_angle_, _init_dose_, _init_a_param_, _init_b_param_, _init_L_alpha_, _init_M_alpha_, _init_rn_on_off_, _init_rn_min_, _init_rn_max_, _init_br_method_, _init_br_family_, _init_br_wavelet_, _init_br_level_, _init_br_max_iterations_, _init_br_order_, _init_br_a0_, _init_br_a1_, _init_br_a2_, _init_br_a3_, _init_br_a4_,  _init_br_lx1_,  _init_br_ly1_,  _init_br_lx2_,  _init_br_ly2_):
            # Load Data
        self.init_file_loc = _init_file_loc_
        self.init_name = _init_name_
        self.init_material = _init_material_
        self.init_crystal_orientation = _init_crystal_orientation_
        self.init_axial_channel = _init_axial_channel_
        self.init_angle = _init_angle_
        self.init_dose = _init_dose_
            # Calibration
        self.init_a_param = _init_a_param_
        self.init_b_param = _init_b_param_
        self.init_L_alpha = _init_L_alpha_
        self.init_M_alpha = _init_M_alpha_
            # Ramove noise 
        self.init_rn_on_off = _init_rn_on_off_
        self.init_rn_min = _init_rn_min_
        self.init_rn_max = _init_rn_max_
            # Background removal
        self.init_br_method = _init_br_method_
        self.init_br_family  = _init_br_family_
        self.init_br_wavelet  = _init_br_wavelet_
        self.init_br_level  = _init_br_level_
        self.init_br_max_iterations  = _init_br_max_iterations_
        self.init_br_order  = _init_br_order_
        self.init_br_a0  = _init_br_a0_
        self.init_br_a1  = _init_br_a1_
        self.init_br_a2  = _init_br_a2_
        self.init_br_a3  = _init_br_a3_
        self.init_br_a4  = _init_br_a4_
        self.init_br_lx1  = _init_br_lx1_
        self.init_br_ly1  = _init_br_ly1_
        self.init_br_lx2  = _init_br_lx2_
        self.init_br_ly2  = _init_br_ly2_

    def save_config(self, file_name): #, init_parameters)
        with open(file_name, "w") as config_file:
            json.dump(self.__dict__, config_file)
        messagebox.showinfo('Info', 'Done!')

    def load_config(file_name):
        with open(file_name, "r") as config_file:
            parameters = json.load(config_file)
            return InitialParameters(parameters.get("init_file_loc"), parameters.get("init_name"), parameters.get("init_material"), parameters.get("init_crystal_orientation"), parameters.get("init_axial_channel"), parameters.get("init_angle"), parameters.get("init_dose"), parameters.get("init_a_param"), parameters.get("init_b_param"), parameters.get("init_L_alpha"), parameters.get("init_M_alpha"), parameters.get("init_rn_on_off"), parameters.get("init_rn_min"), parameters.get("init_rn_max"), parameters.get("init_br_method"), parameters.get("init_br_family"), parameters.get("init_br_wavelet"), parameters.get("init_br_level"), parameters.get("init_br_max_iterations"), parameters.get("init_br_order"), parameters.get("init_br_a0"), parameters.get("init_br_a1"), parameters.get("init_br_a2"), parameters.get("init_br_a3"), parameters.get("init_br_a4"), parameters.get("init_br_lx1"), parameters.get("init_br_ly1"), parameters.get("init_br_lx2"), parameters.get("init_br_ly2"))


#def save_config_parameters(self): 
#    initial_parameters = InitialParameters("Noname", "Random", "10", "4470", "1040", 0, "250", "2500", 0, 0, 0, "8", "150", "15", "0.2", "1.6e-1", "-1e-3", "3.24e-6", "-3.5e-9",  "1525",  "1525",  "1780",  "1780")
#    initial_parameters.save_config("config.json")



