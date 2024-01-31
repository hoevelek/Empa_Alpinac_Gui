# from future import annotations
from __future__ import annotations

import tkinter as tk
from tkinter import ttk
from tkinter import messagebox
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk
from tkinter import Entry, Frame, Label, Tk, filedialog
from pathlib import Path
import json
import os
import numpy as np

# dependencies
from dependencies_for_simple_user_workflow import return_defaults, do_mass_calibration, show_data_raw_RT_plot, show_data_raw_mass_plot, DefaultObjects
from data_analysis_utils import load_data

class MyWidgetApp:

    def get_defaults(self):
        # try to get user-defined defaults from file
        try:
            with open("user_def_defaults.json", "r") as f:
                user_def_defaults = json.load(f)
        except FileNotFoundError:
            # get defaults from return_defaults()
            user_def_defaults = return_defaults()
            #user_def_defaults = user_def_defaults.copy()
            return user_def_defaults

    def set_defaults_as_attributes(self, reset=False):
        if reset == True:
            defaults = return_defaults()
            with open("user_defined_defaults.json", "w") as f:
                json.dump(defaults, f)
        else:
            defaults = self.get_defaults()
        self.defaults_raw_data = DefaultObjects(defaults["defaults_raw_data"])
        self.defaults_mass_cal = DefaultObjects(defaults["defaults_mass_cal"])
        self.defaults_extraction = DefaultObjects(defaults["defaults_extraction"])
        self.defaults_alpinac = DefaultObjects(defaults["defaults_alpinac"])
        self.defaults_matching = DefaultObjects(defaults["defaults_matching"])
        self.user_paths = DefaultObjects(defaults["user_paths"])

    def get_defaults_as_dict_from_attributes(self):
        defaults = {}
        # do backwards conversion from attributes to dict
        defaults["defaults_raw_data"] = DefaultObjects.convert_back_to_dict(self.defaults_raw_data)
        defaults["defaults_mass_cal"] = DefaultObjects.convert_back_to_dict(self.defaults_mass_cal)
        defaults["defaults_extraction"] = DefaultObjects.convert_back_to_dict(self.defaults_extraction)
        defaults["defaults_alpinac"] = DefaultObjects.convert_back_to_dict(self.defaults_alpinac)
        defaults["defaults_matching"] = DefaultObjects.convert_back_to_dict(self.defaults_matching)
        defaults["user_paths"] = DefaultObjects.convert_back_to_dict(self.user_paths)
        return defaults



    def __init__(self, root):

        #self.user_defined_defaults = {}
        #self.user_defined_paths = {}

        # initialize
        # set file_var
        self.file_var_raw = Path()

        self.data_EI = None
        self.data_CI = None
        self.data_EI_reader = None
        self.data_CI_reader = None
        self.ionisation_dict = {}

        #mass cal
        self.path_mass_cal_EI = None
        self.path_mass_cal_CI = None
        self.file_var_mass_cal_file = None
        self.file_var_mass_cal_mode = None

        #extract
        self.extraction_mode = tk.StringVar()

        # get defaults
        self.set_defaults_as_attributes()

        # set message var to ""
        self.message_raw_data_var = tk.StringVar()
        self.message_mass_cal_var = tk.StringVar()
        self.message_extract_var = tk.StringVar()
        self.message_alpinac_var = tk.StringVar()
        self.message_alpinac_input_var = tk.StringVar()
        self.message_matching_var = tk.StringVar()
        

        # Initialize root window
        self.root = root
        self.root.title("Non-target Analysis App")

        # Create main frame with three horizontal sections
        self.main_frame = ttk.Frame(root)
        self.main_frame.pack(expand=1, fill="both")

        # Upper left: Settings menu
        self.create_settings_menu()

        # Upper middle: Retention Time Plot
        self.rt_plot_frame = ttk.Frame(self.main_frame)
        self.rt_plot_frame.grid(row=2, column=1, sticky="nsew")
        self.plot_rt()

        # Lower middle: Mass Plot
        self.mass_plot_frame = ttk.Frame(self.main_frame)
        self.mass_plot_frame.grid(row=3, column=1, sticky="nsew")
        self.plot_mass(file_path="")

        # Workflow section
        self.create_workflow_section()

        # Set row and column weights for resizing
        self.main_frame.rowconfigure(1, weight=1)
        self.main_frame.rowconfigure(2, weight=1)
        self.main_frame.columnconfigure(1, weight=1)


    # SETTINGS

    def create_settings_menu(self):
        # Settings menu
        self.menu_bar = tk.Menu(self.root)
        self.root.config(menu=self.menu_bar)

        settings_menu = tk.Menu(self.menu_bar, tearoff=0)
        self.menu_bar.add_cascade(label="Settings", menu=settings_menu)

        # Defaults submenu
        defaults_menu = tk.Menu(settings_menu, tearoff=0)
        settings_menu.add_cascade(label="Defaults", menu=defaults_menu)
        defaults_menu.add_command(label="Paths", command=self.open_paths_window)
        defaults_menu.add_command(label="Calibrate", command=self.open_calibrate_window)
        defaults_menu.add_command(label="Extract", command=self.open_extract_window)
        defaults_menu.add_command(label="Alpinac", command=self.open_alpinac_window)
        defaults_menu.add_command(label="Matching", command=self.open_matching_window)

        # Reset Defaults option
        settings_menu.add_command(label="Reset Defaults", command=self.reset_defaults)

        # Close submenu
        settings_menu.add_command(label="Close", command=self.root.destroy)


    def open_defaults_window(self, title):
        defaults_window = tk.Toplevel(self.root)
        defaults_window.title(f"{title} Defaults Configuration")
        title_dict_name_conversion = {
            "Paths": self.user_paths,
            "Raw Data": self.defaults_raw_data,
            "Calibrate": self.defaults_mass_cal,
            "Extract": self.defaults_extraction,
            "Alpinac": self.defaults_alpinac,
            "Matching": self.defaults_matching
        }

        # Add Entry widgets for each default value
        entries = {}
        for i, attribute in enumerate(title_dict_name_conversion[title].keys()):
            # if not class attribute neither starting with __
            label = tk.Label(defaults_window, text=f"{attribute}:")
            label.grid(row=i, column=0, padx=5, pady=5)
            entry = tk.Entry(defaults_window, width=100)
            entry.insert(0, getattr(title_dict_name_conversion[title], attribute))
            entry.grid(row=i, column=1, padx=5, pady=5)
            entries[attribute] = entry
            # set attribute to entry
            setattr(title_dict_name_conversion[title], attribute, entry.get())

        # Save button
        save_button = tk.Button(defaults_window, text="Save and Close",
                                command=lambda: self.save_new_defaults(defaults_window))#defaults_window, entries, get_defaults_as_dict_from_attributes()[title_dict_name_conversion[title]]))
        save_button.grid(row=len(entries), column=0, columnspan=1, pady=10)


        # for i, (key, value) in enumerate(dir(title_dict_name_conversion[title].attributes()):
        #     label = tk.Label(defaults_window, text=f"{key}:")
        #     label.grid(row=i, column=0, padx=5, pady=5)

        #     entry = tk.Entry(defaults_window, width=100)
        #     entry.insert(0, value)
        #     entry.grid(row=i, column=1, padx=5, pady=5)
        #     entries[key] = entry

        # # Save button
        # save_button = tk.Button(defaults_window, text="Save and Close",
        #                         command=lambda: self.save_new_defaults(defaults_window, entries, title_dict_name_conversion[title]))
        # save_button.grid(row=len(self.user_defined_defaults[title_dict_name_conversion[title]]), column=0, columnspan=1, pady=10)

    #def save_new_defaults(self, defaults_window, entries, key_title):
    #    for key, entry in entries.items():
    #        self.user_defined_defaults[key_title][key] = entry.get()
    #    # Save user-defined defaults to file
    #    with open("user_defined_defaults.json", "w") as f:
    #        json.dump(self.user_defined_defaults, f)

    def save_new_defaults(self, defaults_window):#, defaults_window, entries, key_title):
        new_defaults = self.get_defaults_as_dict_from_attributes()
        # update new_defaults with user-defined defaults
        # Save user-defined defaults to file
        with open("user_defined_defaults.json", "w") as f:
            json.dump(new_defaults, f)
        # close window
        defaults_window.destroy()


    def open_paths_window(self):
        self.open_defaults_window("Paths")

    def open_raw_data_window(self):
        self.open_defaults_window("Raw Data")

    def open_calibrate_window(self):
        self.open_defaults_window("Calibrate")

    def open_extract_window(self):
        self.open_defaults_window("Extract")

    def open_alpinac_window(self):
        self.open_defaults_window("Alpinac")

    def open_matching_window(self):
        self.open_defaults_window("Matching")

    def reset_defaults(self):
        self.set_defaults_as_attributes(reset=True) #reset defaults (also in user_defined_defaults.json)
        # Reload the windows to reflect the changes
        self.open_paths_window()
        self.open_calibrate_window()
        self.open_extract_window()
        self.open_alpinac_window()
        self.open_matching_window()

    # WORKFLOW

    def create_workflow_section(self):
        # Workflow section
        self.workflow_frame = ttk.Frame(self.main_frame)
        self.workflow_frame.grid(row=0, column=1, sticky="nsew")

        # Raw Data Section
        raw_data_section = ttk.Frame(self.workflow_frame, borderwidth=2, relief="groove")
        raw_data_section.grid(row=0, column=0, padx=5, pady=5)
        self.add_workflow_section(raw_data_section, "Raw Data", "do_import")

        # Mass Calibration Section
        mass_calibration_section = ttk.Frame(self.workflow_frame, borderwidth=2, relief="groove")
        mass_calibration_section.grid(row=0, column=1, padx=5, pady=5)
        self.add_workflow_section(mass_calibration_section, "Mass Calibration", "do_calibration")

        # Extraction Section
        extraction_section = ttk.Frame(self.workflow_frame, borderwidth=2, relief="groove")
        extraction_section.grid(row=0, column=2, padx=5, pady=5)
        self.add_workflow_section(extraction_section, "Extraction", "do_extraction")

        # Identification Alpinac Section
        alpinac_section = ttk.Frame(self.workflow_frame, borderwidth=2, relief="groove")
        alpinac_section.grid(row=0, column=3, padx=5, pady=5)
        self.add_workflow_section(alpinac_section, "Alpinac", "do_alpinac")

        # Import Alpinac Input Section
        #alpinac_input_section = ttk.Frame(self.workflow_frame, borderwidth=2, relief="groove")
        #alpinac_input_section.grid(row=0, column=4, padx=5, pady=5)
        #self.add_workflow_section(alpinac_input_section, "Alpinac Input", "do_alpinac_input")

        # Identification Alpinac Section
        matching_section = ttk.Frame(self.workflow_frame, borderwidth=2, relief="groove")
        matching_section.grid(row=0, column=4, padx=5, pady=5)
        self.add_workflow_section(matching_section, "Matching", "do_matching")


    def browse_file(self):
        file_path = Path(os.path.normpath(filedialog.askopenfilename()))  # Open file dialog
        if file_path.exists():
            # Update message box with file name
            self.message_raw_data_var.set("File Selected: " + str(file_path))
            # Pass the selected file path to plot methods
            #self.plot_rt(file_path)
            #self.plot_mass(file_path)
            self.file_var_raw = file_path
            ei_tuple = tuple(self.defaults_raw_data.EI_segments)
            # remove all entries from EI_segments which are not integers
            ei_tuple = tuple([int(i) for i in ei_tuple if str(i).isnumeric()])
            ci_tuple = tuple(self.defaults_raw_data.CI_segments)
            # remove all entries from CI_segments which are not integers
            ci_tuple = tuple([int(i) for i in ci_tuple if str(i).isnumeric()])
            self.ionisation_dict = {"EI": ei_tuple, "CI": ci_tuple}

            # load data
            if file_path.exists():
                self.message_raw_data_var.set("Valid File selected: {} ".format(str(file_path)))
            else:
                self.message_raw_data_var.set("File {} does can not be loaded".format(str(file_path)))
        else:
            self.message_raw_data_var.set("No Raw Data File Selected")

    def browse_mass_cal_file(self):
        file_path = Path(os.path.normpath(filedialog.askopenfilename()))
        if file_path:
            self.message_mass_cal_var.set("File Selected: " + str(file_path))
            self.file_var_mass_cal_file.set(file_path)
        else:
            self.message_mass_cal_var.set("No Mass Cal File Selected")

    def browse_alpinac_input(self):
        file_path = Path(os.path.normpath(filedialog.askopenfilename()))  # Open file dialog
        if file_path:
            # Update message box with file name
            self.message_alpinac_input_var.set("File Selected: " + str(file_path))
            # Pass the selected file path to plot methods
            # self.plot_rt(file_path)
            # self.plot_mass(file_path)
            self.file_alpinac_input.set(file_path)
            # return file_path
        else:
            self.message_alpinac_input_var.set("No Alpinac Input File Selected")


    def add_workflow_section(self, frame, title, command1, command2=None):
        # Title label
        title_label = tk.Label(frame, text=title, font=("Helvetica", 10, "bold"))
        title_label.grid(row=0, column=0, pady=5)

        # Update Plot button
        if title == "Raw Data":
            self.message_raw_data_var.set("no file selected")
            self.create_message_box_raw_data(frame)

            # Select Files
            # Dropdown for file selection
            #file_dropdown = ttk.Combobox(frame, textvariable=file_var, state="readonly")
            #file_dropdown.grid(row=1, column=1, padx=5, pady=5)
            #file_dropdown["values"] = ["EI", "CI"]  # Example values
            #file_dropdown.current(0)  # Default selection



            # add field called "EI segements" with default value (1,2)
            EI_segments_label = tk.Label(frame, text="EI Segments:")
            EI_segments_label.grid(row=2, column=0, padx=5, pady=5)
            EI_segments_entry = tk.Entry(frame)
            EI_segments_entry.insert(0, self.defaults_raw_data.EI_segments)
            EI_segments_entry.grid(row=2, column=1, padx=5, pady=5)
            self.defaults_raw_data.EI_segments = EI_segments_entry.get()
            
            # same for CI segments
            CI_segments_label = tk.Label(frame, text="CI Segments:")
            CI_segments_label.grid(row=3, column=0, padx=5, pady=5)
            CI_segments_entry = tk.Entry(frame)
            CI_segments_entry.insert(0, self.defaults_raw_data.CI_segments)
            CI_segments_entry.grid(row=3, column=1, padx=5, pady=5)
            self.defaults_raw_data.CI_segments = CI_segments_entry.get()

            ei_tuple = tuple(self.defaults_raw_data.EI_segments)
            # remove all entries from EI_segments which are not integers
            ei_tuple = tuple([int(i) for i in ei_tuple if str(i).isnumeric()])
            ci_tuple = tuple(self.defaults_raw_data.CI_segments)
            # remove all entries from CI_segments which are not integers
            ci_tuple = tuple([int(i) for i in ci_tuple if str(i).isnumeric()])
            self.ionisation_dict = {"EI": ei_tuple, "CI": ci_tuple}
            print("ion_dict is{}:".format(self.ionisation_dict))
            # Button to browse for file
            browse_button = ttk.Button(frame, text="Browse Raw Data", command=self.browse_file)
            browse_button.grid(row=1, column=0, padx=5, pady=5)

            update_plot_button = tk.Button(frame, text="Update Plots", command=self.update_plots)
            update_plot_button.grid(row=1, column=1, columnspan=1, pady=5)

        elif title == "Mass Calibration":
            self.message_mass_cal_var.set("no file selected (required for spiked mode only)")
            self.create_message_box_mass_cal(frame)

            mass_cal_file = tk.StringVar()
            file_dropdown_cal_file = ttk.Combobox(frame, textvariable=mass_cal_file, state="readonly")
            file_dropdown_cal_file.grid(row=2, column=0, padx=5, pady=5)
            file_dropdown_cal_file["values"] = ["EI", "CI"]  # Example values
            # Set default selection
            file_dropdown_cal_file.current(0)
            self.file_var_mass_cal_file = mass_cal_file

            mass_cal_mode = tk.StringVar()
            file_dropdown_cal_mode = ttk.Combobox(frame, textvariable=mass_cal_mode, state="readonly")
            file_dropdown_cal_mode.grid(row=2, column=1, padx=5, pady=5)
            file_dropdown_cal_mode["values"] = ["continous", "spiked"]  # Example values
            # Set default selection
            file_dropdown_cal_mode.current(0)
            if self.defaults_mass_cal.mass_cal_mode == "continous":
                file_dropdown_cal_mode.current(0)
            elif self.defaults_mass_cal.mass_cal_mode == "spiked":
                file_dropdown_cal_mode.current(1)
            self.file_var_mass_cal_mode = mass_cal_mode

            # Button to browse for file
            cal_file_browse_button = ttk.Button(frame, text="Mass Cal. File", command=self.browse_mass_cal_file)
            cal_file_browse_button.grid(row=1, column=0, padx=5, pady=5)

            calibrate_button = tk.Button(frame, text="Calibrate", command=lambda: self.show_message_mass_cal(command1))
            calibrate_button.grid(row=1, column=1, pady=5)

            # add field called "calibrant injected: " with default value "PFTBA"
            calibrant_injected_label = tk.Label(frame, text="Calibrant Injected:")
            calibrant_injected_label.grid(row=3, column=0, padx=5, pady=5)
            calibrant_injected_entry = tk.Entry(frame)
            calibrant_injected_entry.insert(0, self.defaults_mass_cal.calibrant_injected)
            calibrant_injected_entry.grid(row=3, column=1, padx=5, pady=5)
            # save user-defined calibrant_injected to variable
            # self.user_defined_defaults["defaults_mass_cal"]["calibrant_injected"] = calibrant_injected_entry.get()

        elif title == "Extraction":
            self.message_extract_var.set("no raw data file selected")
            self.create_message_box_extract(frame)

            mode_var_extr = tk.StringVar()
            # Select Files
            # Dropdown for file selection
            file_dropdown_extr = ttk.Combobox(frame, textvariable=mode_var_extr, state="readonly")
            file_dropdown_extr.grid(row=1, column=1, padx=5, pady=5)
            file_dropdown_extr["values"] = ["EI", "EI CI", "CI"]  # Example values
            file_dropdown_extr.current(0)  # Default selection
            self.extraction_mode = mode_var_extr

            # Button to browse for file
            #browse_button = ttk.Button(frame, text="Browse Raw Data", command=self.browse_file)
            #browse_button.grid(row=1, column=0, padx=5, pady=5)

            #update_plot_button = tk.Button(frame, text="Update Plots", command=self.update_plots(self.file_var_raw))
            #update_plot_button.grid(row=1, column=1, columnspan=1, pady=5)

            # add field called "EI segements" with default value (1,2)
            RT_start_label = tk.Label(frame, text="RT start:")
            RT_start_label.grid(row=2, column=0, padx=5, pady=5)
            RT_start_entry = tk.Entry(frame)
            RT_start_entry.insert(0, self.defaults_extraction.rt_start_extract)
            RT_start_entry.grid(row=2, column=1, padx=5, pady=5)
            # same for CI segments
            RT_end_label = tk.Label(frame, text="RT end:")
            RT_end_label.grid(row=3, column=0, padx=5, pady=5)
            RT_end_entry = tk.Entry(frame)
            RT_end_entry.insert(0, self.defaults_extraction.rt_stop_extract)
            RT_end_entry.grid(row=3, column=1, padx=5, pady=5)

            extract_button = tk.Button(frame, text="Extract", command=lambda: self.show_message_extract(command1))
            extract_button.grid(row=1, column=0, pady=5)
        elif title == "Alpinac":
            self.message_alpinac_input_var.set("no alpinac input file selected")
            self.create_message_box_alpinac_input(frame)

            self.message_alpinac_var.set("alpinac did not run yet")
            self.create_message_box_alpinac(frame)

            # Button to browse for file
            browse_button = ttk.Button(frame, text="Browse Alpinac Input", command=self.browse_alpinac_input)
            browse_button.grid(row=1, column=0, padx=5, pady=5)

            alpinac_button = tk.Button(frame, text="Alpinac", command=lambda: self.show_message_alpinac(command1))
            alpinac_button.grid(row=1, column=1, pady=5)

            # add field called "Target Elements" with default value self.user_defined_defaults["defaults_alpinac"]["target_elements"]
            target_elements_label = tk.Label(frame, text="Target Elements:")
            target_elements_label.grid(row=2, column=0, padx=5, pady=5)
            target_elements_entry = tk.Entry(frame)
            target_elements_entry.insert(0, self.defaults_alpinac.target_elements)
            target_elements_entry.grid(row=2, column=1, padx=5, pady=5)

            # matching_button = tk.Button(frame, text="Matching", command=lambda: self.show_message(command2))
            # matching_button.grid(row=2, column=0, pady=5)
        elif title == "Matching":
            self.message_matching_var.set("not implemented yet")
            self.create_message_box_matching(frame)
            matching_button = tk.Button(frame, text="Matching", command=lambda: self.show_message_matching(command1))
            matching_button.grid(row=1, column=0, pady=5)


    def show_message_raw_data(self, command):
        # Example implementation to update message box
        result = self.execute_workflow_command(command)
        self.message_raw_data_var.set(result)
        # update message box with file_name

    def show_message_mass_cal(self, command):
        # Example implementation to update message box
        result = self.execute_workflow_command(command)
        self.message_mass_cal_var.set(result)
        # update message box with file_name

    def show_message_extract(self, command):
        # Example implementation to update message box
        result = self.execute_workflow_command(command)
        self.message_extract_var.set(result)
        # update message box with file_name

    def show_message_alpinac_input(self, command):
        # Example implementation to update message box
        result = self.execute_workflow_command(command)
        self.message_alpinac_input_var.set(result)
        # update message box with file_name

    def show_message_alpinac(self, command):
        # Example implementation to update message box
        result = self.execute_workflow_command(command)
        self.message_alpinac_var.set(result)
        # update message box with file_name

    def show_message_matching(self, command):
        # Example implementation to update message box
        result = self.execute_workflow_command(command)
        self.message_matching_var.set(result)
        # update message box with file_name

    def create_message_box_raw_data(self, frame):
        # Message box
        # Example implementation to show messages
        message_box = tk.Entry(frame, textvariable=self.message_raw_data_var, state='readonly', width=38)
        message_box.grid(row=4, column=0, columnspan=2, pady=5)

    def create_message_box_mass_cal(self, frame):
        # Message box
        # Example implementation to show messages
        message_box = tk.Entry(frame, textvariable=self.message_mass_cal_var, state='readonly', width=38)
        message_box.grid(row=4, column=0, columnspan=2, pady=5)

    def create_message_box_extract(self, frame):
        # Message box
        # Example implementation to show messages
        message_box = tk.Entry(frame, textvariable=self.message_extract_var, state='readonly', width=38)
        message_box.grid(row=4, column=0, columnspan=2, pady=5)

    def create_message_box_alpinac(self, frame):
        # Message box
        # Example implementation to show messages
        message_box = tk.Entry(frame, textvariable=self.message_alpinac_var, state='readonly', width=38)
        message_box.grid(row=5, column=0, columnspan=2, pady=5)

    def create_message_box_alpinac_input(self, frame):
        # Message box
        # Example implementation to show messages
        message_box = tk.Entry(frame, textvariable=self.message_alpinac_input_var, state='readonly', width=38)
        message_box.grid(row=4, column=0, columnspan=2, pady=5)

    def create_message_box_matching(self, frame):
        # Message box
        # Example implementation to show messages
        message_box = tk.Entry(frame, textvariable=self.message_matching_var, state='readonly', width=38)
        message_box.grid(row=3, column=0, columnspan=2, pady=5)



    def execute_workflow_command(self, command):
        # Example implementation of workflow commands
        if command == "do_import":
            # Your logic to import raw data
            return "Import successful!"
        elif command == "do_calibration":
            return "Calibration successful!"
        elif command == "do_extraction":
            return "Extraction successful!"
        elif command == "do_alpinac":
            return "Alpinac successful!"
        elif command == "do_matching":
            return "Matching successful!"


    def update_plots(self):
        self.message_raw_data_var.set("Load file {}".format(self.file_var_raw))
        path_var_raw = Path(os.path.normpath(self.file_var_raw))
        ionisation_dict_h = self.ionisation_dict

        # if self.file_var_raw is str replace \ with / and check if file exists
        if isinstance(self.file_var_raw, str):
            self.file_var_raw = Path(self.file_var_raw.replace("\\", "/"))
            if self.file_var_raw.exists():
                self.message_raw_data_var.set("Valid File h5 selected: {} ".format(str(self.file_var_raw)))
        if not Path(os.path.normpath(self.file_var_raw)).exists():
            # print in message box that no file is selected
            self.message_raw_data_var.set("Loading error for file {}".format(self.file_var_raw))
            pass
        else:
            # Example implementation for "Update Plots" button
            print(ionisation_dict_h)
            print(self.file_var_raw)
            print("type of file_var_raw: {}".format(type(self.file_var_raw)))
            print("type of ionisation_dict: {}".format(type(ionisation_dict_h)))
           
            if "EI" in ionisation_dict_h.keys() and path_var_raw.exists():
                print("EI in ionisation_dict_h.keys() and path_var_raw.exists()")
                self.message_raw_data_var.set("Try loading File (EI) {}".format(path_var_raw))
                #path_var_raw = Path(r"C:\Users\kaho\Desktop\data\data_Empa\Campaign202303\230317.2157.tank.5.h5")
                #ionisation_dict_h = {"EI": (3,4), "CI": (1,2)}
                sucess_flag_EI = False
                try:
                    data_EI, reader_EI = load_data(file_path=path_var_raw, ionisation_dict=ionisation_dict_h, ionisation_type="EI")
                    sucess_flag_EI = True
                except:
                    self.message_raw_data_var.set("EI spectra not loaded of File {} ".format(path_var_raw))
                if sucess_flag_EI:
                    self.data_EI = data_EI
                    self.data_EI_reader = reader_EI
                    self.message_raw_data_var.set("Loaded File (EI) {}".format(path_var_raw))

            if "CI" in ionisation_dict_h.keys() and path_var_raw.exists():
                print("CI in ionisation_dict_h.keys() and path_var_raw.exists()")
                sucess_flag_CI = False
                self.message_raw_data_var.set("Try loading File (CI) {}".format(path_var_raw))
                try:
                    data_CI, reader_CI = load_data(path_var_raw, ionisation_dict_h, "CI")
                    #path_var_raw = Path(r"C:\Users\kaho\Desktop\data\data_Empa\Campaign202303\230317.2157.tank.5.h5")
                    sucess_flag_CI = True
                except:
                    self.message_raw_data_var.set("CI spectra not loaded of File {} ".format(path_var_raw))
                if sucess_flag_CI:
                    ionisation_dict_h = {"EI": (3,4), "CI": (1,2)}
                    self.data_CI = data_CI
                    self.data_CI_reader = reader_CI
                    self.message_raw_data_var.set("Loaded File (CI) {}".format(path_var_raw))
            
        
        # if file exists, update plots
        if sucess_flag_CI or sucess_flag_EI:
            self.plot_rt()
        #    self.plot_mass(file_path)

    def plot_rt(self):

        # set value to user_defined_defaults["defaults_raw_data"]["mass_start"]

        # Entry fields for RT start and RT end
        mass_start_label = tk.Label(self.rt_plot_frame, text="Mass Start:")
        mass_start_label.grid(row=1, column=0, padx=5, pady=5)
        mass_start_entry = tk.Entry(self.rt_plot_frame)
        mass_start_entry.insert(0, self.defaults_raw_data.mass_ranges_start)
        mass_start_entry.grid(row=1, column=1, padx=5, pady=5)
        
        # update user_defined_defaults["defaults_raw_data"]["mass_start"] with value from entry field
        self.defaults_raw_data.mass_ranges_start = mass_start_entry.get()

        # set value to user_defined_defaults["defaults_raw_data"]["mass_end"]
        mass_end_label = tk.Label(self.rt_plot_frame, text="Mass End:")
        mass_end_label.grid(row=1, column=2, padx=5, pady=5)
        mass_end_entry = tk.Entry(self.rt_plot_frame)
        mass_end_entry.insert(0, self.defaults_raw_data.mass_ranges_stop)
        mass_end_entry.grid(row=1, column=3, padx=5, pady=5)

        # get value from user entry field



        # update user_defined_defaults["defaults_raw_data"]["mass_start"] and user_defined_defaults["defaults_raw_data"]["mass_end"] with values from entry fields
        self.defaults_raw_data.mass_ranges_stop = mass_end_entry.get()

        # Save button for RT range
        mass_save_button = tk.Button(self.rt_plot_frame, text="Apply Mass Range", command=lambda: self.apply_mass_range(ax, mass_start_entry.get(), mass_end_entry.get()))
        mass_save_button.grid(row=1, column=4, columnspan=1, padx= 5, pady=5)

        
        # get file_path from user entry field
        fig, ax = plt.subplots(figsize=(5, 3), tight_layout=True)
        # test if filepath exist and ends with h5
        #if isinstance(file_path, str) and Path(file_path).exists() and file_path.endswith(".h5"):
        #    ionisation_dict = {"EI": self.user_defined_defaults["defaults_raw_data"]["EI_segments"], "CI": self.user_defined_defaults["defaults_raw_data"]["CI_segments"]}
            # ionisation type: if EI_segments is empty, use CI, if CI_segments is empty, use EI else loop through both
        # if self data Ei is a 2D np array and not empty:
        # get user defined defaults as attributes
        defaults_for_plot = self.get_defaults_as_dict_from_attributes()
        flag_EI_plot = False
        if isinstance(self.data_EI, np.ndarray) and self.data_EI.size > 0:
            print("mass_ranges_start are: {}".format(self.defaults_raw_data.mass_ranges_start))
            print("mass_ranges_stop are: {}".format(self.defaults_raw_data.mass_ranges_stop))
            # convert mass_ranges_start and mass_ranges_stop to str, remove all brackets and split
            mass_ranges_start = str(self.defaults_raw_data.mass_ranges_start).replace("(", "").replace(")", "").split()
            mass_ranges_stop = str(self.defaults_raw_data.mass_ranges_stop).replace("(", "").replace(")", "").split()
            #convert to list
            mass_ranges_start = [int(i) for i in mass_ranges_start]
            mass_ranges_stop = [int(i) for i in mass_ranges_stop]
            self.defaults_raw_data.mass_ranges_start = mass_ranges_start
            self.defaults_raw_data.mass_ranges_stop = mass_ranges_stop
            #print default_for_plots values
            print(self.defaults_raw_data["mass_ranges_start"])
            print(self.defaults_raw_data["mass_ranges_stop"])

            show_data_raw_RT_plot(data=self.data_EI, reader=self.data_EI_reader, defaults=defaults_for_plot, ionisation_dict=self.ionisation_dict, ionisation_type="EI", mass_ranges_start=self.defaults_raw_data.mass_ranges_start, mass_ranges_end=self.defaults_raw_data.mass_ranges_stop, ax=ax)
            #show_data_raw_RT_plot(data=self.data_EI, reader=self.data_EI_reader, defaults=defaults_for_plot, ionisation_dict=self.ionisation_dict, ionisation_type="EI", mass_ranges_start=None, mass_ranges_end=None, ax=ax)
            ax.legend(["EI, mass range: {} - {}".format(self.defaults_raw_data.mass_ranges_start, self.defaults_raw_data.mass_ranges_stop)])            
            flag_EI_plot = True

        flag_CI_plot = False
        if isinstance(self.data_CI, np.ndarray) and self.data_CI.size > 0:
            flag_CI_plot = True
            show_data_raw_RT_plot(data=self.data_CI, reader=self.data_CI_reader, defaults=defaults_for_plot, ionisation_dict=self.ionisation_dict, ionisation_type="CI", mass_ranges_start=self.defaults_raw_data.mass_ranges_start, mass_ranges_end=self.defaults_raw_data.mass_ranges_stop, ax=ax)
            ax.legend(["CI, mass range: {} - {}".format(self.defaults_raw_data.mass_ranges_start, self.defaults_raw_data.mass_ranges_stop)])
            
        if not flag_EI_plot and not flag_CI_plot:
            ax.set_title("Retention Time Plot")
            ax.set_xlabel("Retention Time")
            ax.set_ylabel("Intensity")
            ax.set_title("Retention Time Plot: " + str(self.file_var_raw))

        # # Pack the canvas and toolbar
        canvas = FigureCanvasTkAgg(fig, master=self.rt_plot_frame)
        canvas.draw()
        canvas.get_tk_widget().grid(row=3, column=0, columnspan=5, sticky="nsew")
        #self.mass_plot_frame.columnconfigure(3, weight=1)

        for col in range(4):  # Adjusted: Set weight for all columns
            self.rt_plot_frame.columnconfigure(col, weight=1)

        # Navigation toolbar
        toolbar_frame = tk.Frame(self.rt_plot_frame)
        toolbar = NavigationToolbar2Tk(canvas, toolbar_frame)
        toolbar.update()
        toolbar_frame.grid(row=2, column=0, columnspan=4, padx=5, pady=5, sticky="nsew")

    def apply_mass_range(self, ax, mass_start, mass_end):
        # Example implementation to apply mass range in the plot
        print("Applying Mass Range:", mass_start, mass_end)
        # Your logic to update the plot with the provided mass range
        defaults_for_plot = self.get_defaults_as_dict_from_attributes()
        flag_EI_plot = False
        if isinstance(self.data_EI, np.ndarray) and self.data_EI.size > 0:
            print("mass_ranges_start are: {}".format(self.defaults_raw_data.mass_ranges_start))
            print("mass_ranges_stop are: {}".format(self.defaults_raw_data.mass_ranges_stop))
            # convert mass_ranges_start and mass_ranges_stop to str, remove all brackets and split
            mass_ranges_start = str(self.defaults_raw_data.mass_ranges_start).replace("(", "").replace(")", "").split()
            mass_ranges_stop = str(self.defaults_raw_data.mass_ranges_stop).replace("(", "").replace(")", "").split()
            #convert to list
            mass_ranges_start = [int(i) for i in mass_ranges_start]
            mass_ranges_stop = [int(i) for i in mass_ranges_stop]
            self.defaults_raw_data.mass_ranges_start = mass_ranges_start
            self.defaults_raw_data.mass_ranges_stop = mass_ranges_stop
            #print default_for_plots values
            print(self.defaults_raw_data["mass_ranges_start"])
            print(self.defaults_raw_data["mass_ranges_stop"])

            show_data_raw_RT_plot(data=self.data_EI, reader=self.data_EI_reader, defaults=defaults_for_plot, ionisation_dict=self.ionisation_dict, ionisation_type="EI", mass_ranges_start=self.defaults_raw_data.mass_ranges_start, mass_ranges_end=self.defaults_raw_data.mass_ranges_stop, ax=ax)
            #show_data_raw_RT_plot(data=self.data_EI, reader=self.data_EI_reader, defaults=defaults_for_plot, ionisation_dict=self.ionisation_dict, ionisation_type="EI", mass_ranges_start=None, mass_ranges_end=None, ax=ax)
            # add legend EI + mass range
            ax.legend(["EI, mass range: {} - {}".format(self.defaults_raw_data.mass_ranges_start, self.defaults_raw_data.mass_ranges_stop)])            
            flag_EI_plot = True
        flag_CI_plot = False
        if isinstance(self.data_CI, np.ndarray) and self.data_CI.size > 0:
            flag_CI_plot = True
            show_data_raw_RT_plot(data=self.data_CI, reader=self.data_CI_reader, defaults=defaults_for_plot, ionisation_dict=self.ionisation_dict, ionisation_type="CI", mass_ranges_start=self.defaults_raw_data.mass_ranges_start, mass_ranges_end=self.defaults_raw_data.mass_ranges_stop, ax=ax)
                # add legend EI + mass range
            ax.legend(["CI, mass range: {} - {}".format(self.defaults_raw_data.mass_ranges_start, self.defaults_raw_data.mass_ranges_stop)])            
        if not flag_EI_plot and not flag_CI_plot:
            ax.set_title("Retention Time Plot")
            ax.set_xlabel("Retention Time")
            ax.set_ylabel("Intensity")
            ax.set_title("Retention Time Plot: " + str(self.file_var_raw))




    def plot_mass(self, file_path):
        # Example implementation for Mass Plot
        fig, ax = plt.subplots(figsize=(5, 3), tight_layout=True)
        ax.set_title("Mass Plot")
        ax.set_xlabel("Mass")
        ax.set_ylabel("Intensity")
        ax.set_title("Mass Plot: " + str(file_path))

        # Entry fields for RT start and RT end
        rt_start_label = tk.Label(self.mass_plot_frame, text="RT Start:")
        rt_start_label.grid(row=1, column=0, padx=5, pady=5)
        rt_start_entry = tk.Entry(self.mass_plot_frame)
        rt_start_entry.grid(row=1, column=1, padx=5, pady=5)

        rt_end_label = tk.Label(self.mass_plot_frame, text="RT End:")
        rt_end_label.grid(row=1, column=2, padx=5, pady=5)
        rt_end_entry = tk.Entry(self.mass_plot_frame)
        rt_end_entry.grid(row=1, column=3, padx=5, pady=5)

        # Save button for RT range
        rt_save_button = tk.Button(self.mass_plot_frame, text="Apply RT Range", command=lambda: self.apply_rt_range(ax, rt_start_entry.get(), rt_end_entry.get()))
        rt_save_button.grid(row=1, column=4, columnspan=1, padx= 5, pady=5)

        # # Pack the canvas and toolbar
        canvas = FigureCanvasTkAgg(fig, master=self.mass_plot_frame)
        canvas.draw()
        canvas.get_tk_widget().grid(row=3, column=0, columnspan=5, sticky="nsew")
        #self.mass_plot_frame.columnconfigure(3, weight=1)

        for col in range(4):  # Adjusted: Set weight for all columns
            self.mass_plot_frame.columnconfigure(col, weight=1)

        # Navigation toolbar
        toolbar_frame = tk.Frame(self.mass_plot_frame)
        toolbar = NavigationToolbar2Tk(canvas, toolbar_frame)
        toolbar.update()
        toolbar_frame.grid(row=2, column=0, columnspan=4, padx=5, pady=5, sticky="nsew")
   
    def apply_rt_range(self, ax, rt_start, rt_end):
        # Example implementation to apply RT range in the plot
        print("Applying RT Range:", rt_start, rt_end)
        # Your logic to update the plot with the provided RT range

if __name__ == "__main__":
    root = tk.Tk()
    app = MyWidgetApp(root)
    root.mainloop()

#import os
#path_str = Path(r"C:/Users/kaho/Desktop/data/data_Empa/Campaign202303/230311.1455.tank.9.h5").parent
# path_str = Path(os.path.normpath("C:/Users/kaho/Desktop/data/data_Empa/Campaign202303/230311.1455.tank.9.h5"))
# path_str2 = Path(r str("C:\Users\kaho\Desktop\data\data_Empa\Campaign202303\230316.0650.tank.11.h5"))
path_str = Path(r"C:\Users\kaho\Desktop\data\data_Empa\Campaign202303\230317.2157.tank.5.h5")

returns = load_data(path_str, {"EI": (3,4), "CI": (1,2)}, "EI")
data, reader = returns