import tkinter as tk
from tkinter import ttk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

class MyWidgetApp:
    def __init__(self, root):
        self.root = root
        self.root.title("My Widget App")

        # Create main frame with three horizontal sections
        self.main_frame = ttk.Frame(root)
        self.main_frame.pack(expand=1, fill="both")

        # Upper left: Settings menu
        self.create_settings_menu()

        # Upper middle: Retention Time Plot
        self.rt_plot_frame = ttk.Frame(self.main_frame)
        self.rt_plot_frame.grid(row=0, column=1, sticky="nsew")
        self.plot_rt()

        # Lower middle: Mass Plot
        self.mass_plot_frame = ttk.Frame(self.main_frame)
        self.mass_plot_frame.grid(row=1, column=1, sticky="nsew")
        self.plot_mass()

        # Set row and column weights for resizing
        self.main_frame.rowconfigure(0, weight=1)
        self.main_frame.rowconfigure(1, weight=1)
        self.main_frame.columnconfigure(1, weight=1)

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
        defaults_menu.add_command(label="Identify", command=self.open_identify_window)

        # Close submenu
        settings_menu.add_command(label="Close", command=self.root.destroy)

    def open_paths_window(self):
        # Example implementation for "Paths" window
        paths_window = tk.Toplevel(self.root)
        paths_window.title("Paths Configuration")

        # Add Entry widgets and Dropdowns for each path
        path1_label = tk.Label(paths_window, text="Path 1:")
        path1_label.grid(row=0, column=0, padx=5, pady=5)
        path1_entry = tk.Entry(paths_window)
        path1_entry.grid(row=0, column=1, padx=5, pady=5)

        # Similar entries and dropdowns for other paths

        # Save button
        save_button = tk.Button(paths_window, text="Save and Close", command=paths_window.destroy)
        save_button.grid(row=1, column=0, columnspan=2, pady=10)

    def open_calibrate_window(self):
        # Example implementation for "Calibrate Defaults" window
        calibrate_window = tk.Toplevel(self.root)
        calibrate_window.title("Calibrate Defaults Configuration")

        # Add Entry widgets for each calibrate default
        calibrate_default1_label = tk.Label(calibrate_window, text="Calibrate Default 1:")
        calibrate_default1_label.grid(row=0, column=0, padx=5, pady=5)
        calibrate_default1_entry = tk.Entry(calibrate_window)
        calibrate_default1_entry.grid(row=0, column=1, padx=5, pady=5)

        # Similar entries for other calibrate defaults

        # Save button
        save_button = tk.Button(calibrate_window, text="Save and Close", command=calibrate_window.destroy)
        save_button.grid(row=1, column=0, columnspan=2, pady=10)

    def open_extract_window(self):
        # Example implementation for "Extract Defaults" window
        extract_window = tk.Toplevel(self.root)
        extract_window.title("Extract Defaults Configuration")

        # Add Entry widgets for each extract default
        extract_default1_label = tk.Label(extract_window, text="Extract Default 1:")
        extract_default1_label.grid(row=0, column=0, padx=5, pady=5)
        extract_default1_entry = tk.Entry(extract_window)
        extract_default1_entry.grid(row=0, column=1, padx=5, pady=5)

        # Similar entries for other extract defaults

        # Save button
        save_button = tk.Button(extract_window, text="Save and Close", command=extract_window.destroy)
        save_button.grid(row=1, column=0, columnspan=2, pady=10)

    def open_identify_window(self):
        # Example implementation for "Identify Defaults" window
        identify_window = tk.Toplevel(self.root)
        identify_window.title("Identify Defaults Configuration")

        # Add Entry widgets for each identify default
        identify_default1_label = tk.Label(identify_window, text="Identify Default 1:")
        identify_default1_label.grid(row=0, column=0, padx=5, pady=5)
        identify_default1_entry = tk.Entry(identify_window)
        identify_default1_entry.grid(row=0, column=1, padx=5, pady=5)

        # Similar entries for other identify defaults

        # Save button
        save_button = tk.Button(identify_window, text="Save and Close", command=identify_window.destroy)
        save_button.grid(row=1, column=0, columnspan=2, pady=10)

    def plot_rt(self):
        # Example implementation for Retention Time Plot
        fig, ax = plt.subplots(figsize=(5, 3), tight_layout=True)
        ax.set_title("Retention Time Plot")
        canvas = FigureCanvasTkAgg(fig, master=self.rt_plot_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill="both", expand=True)

    def plot_mass(self):
        # Example implementation for Mass Plot
        fig, ax = plt.subplots(figsize=(5, 3), tight_layout=True)
        ax.set_title("Mass Plot")
        canvas = FigureCanvasTkAgg(fig, master=self.mass_plot_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill="both", expand=True)

if __name__ == "__main__":
    root = tk.Tk()
    app = MyWidgetApp(root)
    root.mainloop()



import tkinter as tk
from tkinter import ttk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

class MyWidgetApp:
    def __init__(self, root):
        self.root = root
        self.root.title("My Widget App")

        # Create main frame with three horizontal sections
        self.main_frame = ttk.Frame(root)
        self.main_frame.pack(expand=1, fill="both")

        # Upper left: Settings menu
        self.create_settings_menu()

        # Upper middle: Retention Time Plot
        self.rt_plot_frame = ttk.Frame(self.main_frame)
        self.rt_plot_frame.grid(row=0, column=1, sticky="nsew")
        self.plot_rt()

        # Lower middle: Mass Plot
        self.mass_plot_frame = ttk.Frame(self.main_frame)
        self.mass_plot_frame.grid(row=1, column=1, sticky="nsew")
        self.plot_mass()

        # Set row and column weights for resizing
        self.main_frame.rowconfigure(0, weight=1)
        self.main_frame.rowconfigure(1, weight=1)
        self.main_frame.columnconfigure(1, weight=1)

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
        defaults_menu.add_command(label="Identify", command=self.open_identify_window)

        # Close submenu
        settings_menu.add_command(label="Close", command=self.root.destroy)

    def open_paths_window(self):
        # Example implementation for "Paths" window
        paths_window = tk.Toplevel(self.root)
        paths_window.title("Paths Configuration")

        # Add Entry widgets and Dropdowns for each path
        path1_label = tk.Label(paths_window, text="Path 1:")
        path1_label.grid(row=0, column=0, padx=5, pady=5)
        path1_entry = tk.Entry(paths_window)
        path1_entry.grid(row=0, column=1, padx=5, pady=5)

        # Similar entries and dropdowns for other paths

        # Save button
        save_button = tk.Button(paths_window, text="Save and Close", command=paths_window.destroy)
        save_button.grid(row=1, column=0, columnspan=2, pady=10)

    def open_calibrate_window(self):
        # Example implementation for "Calibrate Defaults" window
        calibrate_window = tk.Toplevel(self.root)
        calibrate_window.title("Calibrate Defaults Configuration")

        # Add Entry widgets for each calibrate default
        calibrate_default1_label = tk.Label(calibrate_window, text="Calibrate Default 1:")
        calibrate_default1_label.grid(row=0, column=0, padx=5, pady=5)
        calibrate_default1_entry = tk.Entry(calibrate_window)
        calibrate_default1_entry.grid(row=0, column=1, padx=5, pady=5)

        # Similar entries for other calibrate defaults

        # Save button
        save_button = tk.Button(calibrate_window, text="Save and Close", command=calibrate_window.destroy)
        save_button.grid(row=1, column=0, columnspan=2, pady=10)

    def open_extract_window(self):
        # Example implementation for "Extract Defaults" window
        extract_window = tk.Toplevel(self.root)
        extract_window.title("Extract Defaults Configuration")

        # Add Entry widgets for each extract default
        extract_default1_label = tk.Label(extract_window, text="Extract Default 1:")
        extract_default1_label.grid(row=0, column=0, padx=5, pady=5)
        extract_default1_entry = tk.Entry(extract_window)
        extract_default1_entry.grid(row=0, column=1, padx=5, pady=5)

        # Similar entries for other extract defaults

        # Save button
        save_button = tk.Button(extract_window, text="Save and Close", command=extract_window.destroy)
        save_button.grid(row=1, column=0, columnspan=2, pady=10)

    def open_identify_window(self):
        # Example implementation for "Identify Defaults" window
        identify_window = tk.Toplevel(self.root)
        identify_window.title("Identify Defaults Configuration")

        # Add Entry widgets for each identify default
        identify_default1_label = tk.Label(identify_window, text="Identify Default 1:")
        identify_default1_label.grid(row=0, column=0, padx=5, pady=5)
        identify_default1_entry = tk.Entry(identify_window)
        identify_default1_entry.grid(row=0, column=1, padx=5, pady=5)

        # Similar entries for other identify defaults

        # Save button
        save_button = tk.Button(identify_window, text="Save and Close", command=identify_window.destroy)
        save_button.grid(row=1, column=0, columnspan=2, pady=10)

    def plot_rt(self):
        # Example implementation for Retention Time Plot
        fig, ax = plt.subplots(figsize=(5, 3), tight_layout=True)
        ax.set_title("Retention Time Plot")

        # Entry fields for mass start and mass end
        mass_start_label = tk.Label(self.rt_plot_frame, text="Mass Start:")
        mass_start_label.grid(row=0, column=0, padx=5, pady=5)
        self.mass_start_entry = tk.Entry(self.rt_plot_frame)
        self.mass_start_entry.grid(row=0, column=1, padx=5, pady=5)

        mass_end_label = tk.Label(self.rt_plot_frame, text="Mass End:")
        mass_end_label.grid(row=0, column=2, padx=5, pady=5)
        self.mass_end_entry = tk.Entry(self.rt_plot_frame)
        self.mass_end_entry.grid(row=0, column=3, padx=5, pady=5)

        canvas = FigureCanvasTkAgg(fig, master=self.rt_plot_frame)
        canvas.draw()
        canvas.get_tk_widget().grid(row=1, column=0, columnspan=4, sticky="nsew")

    def plot_mass(self):
        # Example implementation for Mass Plot
        fig, ax = plt.subplots(figsize=(5, 3), tight_layout=True)
        ax.set_title("Mass Plot")

        # Entry fields for RT start and RT end
        rt_start_label = tk.Label(self.mass_plot_frame, text="RT Start:")
        rt_start_label.grid(row=0, column=0, padx=5, pady=5)
        self.rt_start_entry = tk.Entry(self.mass_plot_frame)
        self.rt_start_entry.grid(row=0, column=1, padx=5, pady=5)

        rt_end_label = tk.Label(self.mass_plot_frame, text="RT End:")
        rt_end_label.grid(row=0, column=2, padx=5, pady=5)
        self.rt_end_entry = tk.Entry(self.mass_plot_frame)
        self.rt_end_entry.grid(row=0, column=3, padx=5, pady=5)

        canvas = FigureCanvasTkAgg(fig, master=self.mass_plot_frame)
        canvas.draw()
        canvas.get_tk_widget().grid(row=1, column=0, columnspan=4, sticky="nsew")

if __name__ == "__main__":
    root = tk.Tk()
    app = MyWidgetApp(root)
    root.mainloop()



import tkinter as tk
from tkinter import ttk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

class MyWidgetApp:
    def __init__(self, root):
        self.root = root
        self.root.title("My Widget App")

        # Create main frame with three horizontal sections
        self.main_frame = ttk.Frame(root)
        self.main_frame.pack(expand=1, fill="both")

        # Upper left: Settings menu
        self.create_settings_menu()

        # Upper middle: Retention Time Plot
        self.rt_plot_frame = ttk.Frame(self.main_frame)
        self.rt_plot_frame.grid(row=0, column=1, sticky="nsew")
        self.plot_rt()

        # Lower middle: Mass Plot
        self.mass_plot_frame = ttk.Frame(self.main_frame)
        self.mass_plot_frame.grid(row=1, column=1, sticky="nsew")
        self.plot_mass()

        # Set row and column weights for resizing
        self.main_frame.rowconfigure(0, weight=1)
        self.main_frame.rowconfigure(1, weight=1)
        self.main_frame.columnconfigure(1, weight=1)

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
        defaults_menu.add_command(label="Identify", command=self.open_identify_window)

        # Close submenu
        settings_menu.add_command(label="Close", command=self.root.destroy)

    def open_paths_window(self):
        # Example implementation for "Paths" window
        paths_window = tk.Toplevel(self.root)
        paths_window.title("Paths Configuration")

        # Add Entry widgets and Dropdowns for each path
        path1_label = tk.Label(paths_window, text="Path 1:")
        path1_label.grid(row=0, column=0, padx=5, pady=5)
        path1_entry = tk.Entry(paths_window)
        path1_entry.grid(row=0, column=1, padx=5, pady=5)

        # Similar entries and dropdowns for other paths

        # Save button
        save_button = tk.Button(paths_window, text="Save and Close", command=paths_window.destroy)
        save_button.grid(row=1, column=0, columnspan=2, pady=10)

    def open_calibrate_window(self):
        # Example implementation for "Calibrate Defaults" window
        calibrate_window = tk.Toplevel(self.root)
        calibrate_window.title("Calibrate Defaults Configuration")

        # Add Entry widgets for each calibrate default
        calibrate_default1_label = tk.Label(calibrate_window, text="Calibrate Default 1:")
        calibrate_default1_label.grid(row=0, column=0, padx=5, pady=5)
        calibrate_default1_entry = tk.Entry(calibrate_window)
        calibrate_default1_entry.grid(row=0, column=1, padx=5, pady=5)

        # Similar entries for other calibrate defaults

        # Save button
        save_button = tk.Button(calibrate_window, text="Save and Close", command=calibrate_window.destroy)
        save_button.grid(row=1, column=0, columnspan=2, pady=10)

    def open_extract_window(self):
        # Example implementation for "Extract Defaults" window
        extract_window = tk.Toplevel(self.root)
        extract_window.title("Extract Defaults Configuration")

        # Add Entry widgets for each extract default
        extract_default1_label = tk.Label(extract_window, text="Extract Default 1:")
        extract_default1_label.grid(row=0, column=0, padx=5, pady=5)
        extract_default1_entry = tk.Entry(extract_window)
        extract_default1_entry.grid(row=0, column=1, padx=5, pady=5)

        # Similar entries for other extract defaults

        # Save button
        save_button = tk.Button(extract_window, text="Save and Close", command=extract_window.destroy)
        save_button.grid(row=1, column=0, columnspan=2, pady=10)

    def open_identify_window(self):
        # Example implementation for "Identify Defaults" window
        identify_window = tk.Toplevel(self.root)
        identify_window.title("Identify Defaults Configuration")

        # Add Entry widgets for each identify default
        identify_default1_label = tk.Label(identify_window, text="Identify Default 1:")
        identify_default1_label.grid(row=0, column=0, padx=5, pady=5)
        identify_default1_entry = tk.Entry(identify_window)
        identify_default1_entry.grid(row=0, column=1, padx=5, pady=5)

        # Similar entries for other identify defaults

        # Save button
        save_button = tk.Button(identify_window, text="Save and Close", command=identify_window.destroy)
        save_button.grid(row=1, column=0, columnspan=2, pady=10)

    def plot_rt(self):
        # Example implementation for Retention Time Plot
        fig, ax = plt.subplots(figsize=(5, 3), tight_layout=True)
        ax.set_title("Retention Time Plot")

        # Entry fields for mass start and mass end
        mass_start_label = tk.Label(self.rt_plot_frame, text="Mass Start:")
        mass_start_label.grid(row=0, column=0, padx=5, pady=5)
        self.mass_start_entry = tk.Entry(self.rt_plot_frame)
        self.mass_start_entry.grid(row=0, column=1, padx=5, pady=5)

        mass_end_label = tk.Label(self.rt_plot_frame, text="Mass End:")
        mass_end_label.grid(row=0, column=2, padx=5, pady=5)
        self.mass_end_entry = tk.Entry(self.rt_plot_frame)
        self.mass_end_entry.grid(row=0, column=3, padx=5, pady=5)

        canvas = FigureCanvasTkAgg(fig, master=self.rt_plot_frame)
        canvas.draw()
        canvas.get_tk_widget().grid(row=1, column=0, columnspan=4, sticky="nsew")

        # Configure column to expand with the window
        self.rt_plot_frame.columnconfigure(0, weight=1)

    def plot_mass(self):
        # Example implementation for Mass Plot
        fig, ax = plt.subplots(figsize=(5, 3), tight_layout=True)
        ax.set_title("Mass Plot")

        # Entry fields for RT start and RT end
        rt_start_label = tk.Label(self.mass_plot_frame, text="RT Start:")
        rt_start_label.grid(row=0, column=0, padx=5, pady=5)
        self.rt_start_entry = tk.Entry(self.mass_plot_frame)
        self.rt_start_entry.grid(row=0, column=1, padx=5, pady=5)

        rt_end_label = tk.Label(self.mass_plot_frame, text="RT End:")
        rt_end_label.grid(row=0, column=2, padx=5, pady=5)
        self.rt_end_entry = tk.Entry(self.mass_plot_frame)
        self.rt_end_entry.grid(row=0, column=3, padx=5, pady=5)

        canvas = FigureCanvasTkAgg(fig, master=self.mass_plot_frame)
        canvas.draw()
        canvas.get_tk_widget().grid(row=1, column=0, columnspan=4, sticky="nsew")

        # Configure column to expand with the window
        self.mass_plot_frame.columnconfigure(0, weight=1)

if __name__ == "__main__":
    root = tk.Tk()
    app = MyWidgetApp(root)
    root.mainloop()



import tkinter as tk
from tkinter import ttk
from tkinter import messagebox
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

class MyWidgetApp:
    def __init__(self, root):
        self.root = root
        self.root.title("My Widget App")

        # Create main frame with three horizontal sections
        self.main_frame = ttk.Frame(root)
        self.main_frame.pack(expand=1, fill="both")

        # Upper left: Settings menu
        self.create_settings_menu()

        # Workflow section
        self.create_workflow_section()

        # Upper middle: Retention Time Plot
        self.rt_plot_frame = ttk.Frame(self.main_frame)
        self.rt_plot_frame.grid(row=1, column=1, sticky="nsew")
        self.plot_rt()

        # Lower middle: Mass Plot
        self.mass_plot_frame = ttk.Frame(self.main_frame)
        self.mass_plot_frame.grid(row=2, column=1, sticky="nsew")
        self.plot_mass()

        # Set row and column weights for resizing
        self.main_frame.rowconfigure(1, weight=1)
        self.main_frame.rowconfigure(2, weight=1)
        self.main_frame.columnconfigure(1, weight=1)

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
        defaults_menu.add_command(label="Identify", command=self.open_identify_window)

        # Close submenu
        settings_menu.add_command(label="Close", command=self.root.destroy)

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

        # Identification Section
        identification_section = ttk.Frame(self.workflow_frame, borderwidth=2, relief="groove")
        identification_section.grid(row=0, column=3, padx=5, pady=5)
        self.add_workflow_section(identification_section, "Identification", "do_alpinac", "do_matching")

    def add_workflow_section(self, frame, title, command1, command2=None):
        # Title label
        title_label = tk.Label(frame, text=title, font=("Helvetica", 10, "bold"))
        title_label.grid(row=0, column=0, pady=5)

        # Update Plot button
        update_plot_button = tk.Button(frame, text="Update Plot", command=self.update_plot)
        update_plot_button.grid(row=1, column=0, pady=5)

        # Calibrate button
        calibrate_button = tk.Button(frame, text="Calibrate", command=lambda: self.show_message(command1))
        calibrate_button.grid(row=2, column=0, pady=5)

        # Extract button
        extract_button = tk.Button(frame, text="Extract", command=lambda: self.show_message(command1))
        extract_button.grid(row=3, column=0, pady=5)

        # Alpinac button
        alpinac_button = tk.Button(frame, text="Alpinac", command=lambda: self.show_message(command1))
        alpinac_button.grid(row=4, column=0, pady=5)

        # Matching button (only for Identification Section)
        if command2:
            matching_button = tk.Button(frame, text="Matching", command=lambda: self.show_message(command2))
            matching_button.grid(row=5, column=0, pady=5)

        # Message box
        self.create_message_box(frame)

    def create_message_box(self, frame):
        # Message box
        message_var = tk.StringVar()
        message_box = tk.Entry(frame, textvariable=message_var, state='readonly', width=20)
        message_box.grid(row=6, column=0, pady=5)

    def show_message(self, command):
        # Example implementation to update message box
        result = self.execute_workflow_command(command)
        message_var.set(result)

    def execute_workflow_command(self, command):
        # Example implementation of workflow commands
        if command == "do_import":
            return "Raw data imported successfully!"
        elif command == "do_calibration":
            return "Calibration successful!"
        elif command == "do_extraction":
            return "Extraction successful!"
        elif command == "do_alpinac":
            return "Alpinac successful!"
        elif command == "do_matching":
            return "Matching successful!"

    def open_paths_window(self):
        # Example implementation for "Paths" window
        paths_window = tk.Toplevel(self.root)
        paths_window.title("Paths Configuration")

        # Add Entry widgets and Dropdowns for each path
        path1_label = tk.Label(paths_window, text="Path 1:")
        path1_label.grid(row=0, column=0, padx=5, pady=5)
        path1_entry = tk.Entry(paths_window)
        path1_entry.grid(row=0, column=1, padx=5, pady=5)

        # Similar entries and dropdowns for other paths

        # Save button
        save_button = tk.Button(paths_window, text="Save and Close", command=paths_window.destroy)
        save_button.grid(row=1, column=0, columnspan=2, pady=10)

    def open_calibrate_window(self):
        # Example implementation for "Calibrate Defaults" window
        calibrate_window = tk.Toplevel(self.root)
        calibrate_window.title("Calibrate Defaults Configuration")

        # Add Entry widgets for each calibrate default
        calibrate_default1_label = tk.Label(calibrate_window, text="Calibrate Default 1:")
        calibrate_default1_label.grid(row=0, column=0, padx=5, pady=5)
        calibrate_default1_entry = tk.Entry(calibrate_window)
        calibrate_default1_entry.grid(row=0, column=1, padx=5, pady=5)

        # Similar entries for other calibrate defaults

        # Save button
        save_button = tk.Button(calibrate_window, text="Save and Close", command=calibrate_window.destroy)
        save_button.grid(row=1, column=0, columnspan=2, pady=10)

    def open_extract_window(self):
        # Example implementation for "Extract Defaults" window
        extract_window = tk.Toplevel(self.root)
        extract_window.title("Extract Defaults Configuration")

        # Add Entry widgets for each extract default
        extract_default1_label = tk.Label(extract_window, text="Extract Default 1:")
        extract_default1_label.grid(row=0, column=0, padx=5, pady=5)
        extract_default1_entry = tk.Entry(extract_window)
        extract_default1_entry.grid(row=0, column=1, padx=5, pady=5)

        # Similar entries for other extract defaults

        # Save button
        save_button = tk.Button(extract_window, text="Save and Close", command=extract_window.destroy)
        save_button.grid(row=1, column=0, columnspan=2, pady=10)

    def open_identify_window(self):
        # Example implementation for "Identify Defaults" window
        identify_window = tk.Toplevel(self.root)
        identify_window.title("Identify Defaults Configuration")

        # Add Entry widgets for each identify default
        identify_default1_label = tk.Label(identify_window, text="Identify Default 1:")
        identify_default1_label.grid(row=0, column=0, padx=5, pady=5)
        identify_default1_entry = tk.Entry(identify_window)
        identify_default1_entry.grid(row=0, column=1, padx=5, pady=5)

        # Similar entries for other identify defaults

        # Save button
        save_button = tk.Button(identify_window, text="Save and Close", command=identify_window.destroy)
        save_button.grid(row=1, column=0, columnspan=2, pady=10)

    def update_plot(self):
        # Example implementation for "Update Plot" button
        self.plot_rt()
        self.plot_mass()

    def plot_rt(self):
        # Example implementation for Retention Time Plot
        fig, ax = plt.subplots(figsize=(5, 3), tight_layout=True)
        ax.set_title("Retention Time Plot")
        ax.set_xlabel("Retention Time")
        ax.set_ylabel("Intensity")
        canvas = FigureCanvasTkAgg(fig, master=self.rt_plot_frame)
        canvas.draw()
        canvas.get_tk_widget().grid(row=1, column=0, sticky="nsew")
        self.rt_plot_frame.columnconfigure(0, weight=1)

    def plot_mass(self):
        # Example implementation for Mass Plot
        fig, ax = plt.subplots(figsize=(5, 3), tight_layout=True)
        ax.set_title("Mass Plot")
        ax.set_xlabel("Mass")
        ax.set_ylabel("Intensity")
        canvas = FigureCanvasTkAgg(fig, master=self.mass_plot_frame)
        canvas.draw()
        canvas.get_tk_widget().grid(row=1, column=0, sticky="nsew")
        self.mass_plot_frame.columnconfigure(0, weight=1)

if __name__ == "__main__":
    root = tk.Tk()
    app = MyWidgetApp(root)
    root.mainloop()



import tkinter as tk
from tkinter import ttk
from tkinter import messagebox
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

class MyWidgetApp:
    def __init__(self, root):
        self.root = root
        self.root.title("My Widget App")

        # Create main frame with three horizontal sections
        self.main_frame = ttk.Frame(root)
        self.main_frame.pack(expand=1, fill="both")

        # Upper left: Settings menu
        self.create_settings_menu()

        # Workflow section
        self.create_workflow_section()

        # Upper middle: Retention Time Plot
        self.rt_plot_frame = ttk.Frame(self.main_frame)
        self.rt_plot_frame.grid(row=1, column=1, sticky="nsew")
        self.plot_rt()

        # Lower middle: Mass Plot
        self.mass_plot_frame = ttk.Frame(self.main_frame)
        self.mass_plot_frame.grid(row=2, column=1, sticky="nsew")
        self.plot_mass()

        # Set row and column weights for resizing
        self.main_frame.rowconfigure(1, weight=1)
        self.main_frame.rowconfigure(2, weight=1)
        self.main_frame.columnconfigure(1, weight=1)

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
        defaults_menu.add_command(label="Identify", command=self.open_identify_window)

        # Close submenu
        settings_menu.add_command(label="Close", command=self.root.destroy)

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

        # Identification Section
        identification_section = ttk.Frame(self.workflow_frame, borderwidth=2, relief="groove")
        identification_section.grid(row=0, column=3, padx=5, pady=5)
        self.add_workflow_section(identification_section, "Identification", "do_alpinac", "do_matching")

    def add_workflow_section(self, frame, title, command1, command2=None):
        # Title label
        title_label = tk.Label(frame, text=title, font=("Helvetica", 10, "bold"))
        title_label.grid(row=0, column=0, pady=5)

        # Update Plot button
        if title == "Raw Data":
            update_plot_button = tk.Button(frame, text="Update Plots", command=self.update_plots)
            update_plot_button.grid(row=1, column=0, pady=5)
        elif title == "Mass Calibration":
            calibrate_button = tk.Button(frame, text="Calibrate", command=lambda: self.show_message(command1))
            calibrate_button.grid(row=1, column=0, pady=5)
        elif title == "Extraction":
            extract_button = tk.Button(frame, text="Extract", command=lambda: self.show_message(command1))
            extract_button.grid(row=1, column=0, pady=5)
        elif title == "Identification":
            alpinac_button = tk.Button(frame, text="Alpinac", command=lambda: self.show_message(command1))
            alpinac_button.grid(row=1, column=0, pady=5)
            matching_button = tk.Button(frame, text="Matching", command=lambda: self.show_message(command2))
            matching_button.grid(row=2, column=0, pady=5)

        # Message box
        self.create_message_box(frame)

    def create_message_box(self, frame):
        # Message box
        message_var = tk.StringVar()
        message_box = tk.Entry(frame, textvariable=message_var, state='readonly', width=20)
        message_box.grid(row=3, column=0, pady=5)

    def show_message(self, command):
        # Example implementation to update message box
        result = self.execute_workflow_command(command)
        message_var.set(result)

    def execute_workflow_command(self, command):
        # Example implementation of workflow commands
        if command == "do_import":
            return "Raw data imported successfully!"
        elif command == "do_calibration":
            return "Calibration successful!"
        elif command == "do_extraction":
            return "Extraction successful!"
        elif command == "do_alpinac":
            return "Alpinac successful!"
        elif command == "do_matching":
            return "Matching successful!"

    def open_paths_window(self):
        # Example implementation for "Paths" window
        paths_window = tk.Toplevel(self.root)
        paths_window.title("Paths Configuration")

        # Add Entry widgets and Dropdowns for each path
        path1_label = tk.Label(paths_window, text="Path 1:")
        path1_label.grid(row=0, column=0, padx=5, pady=5)
        path1_entry = tk.Entry(paths_window)
        path1_entry.grid(row=0, column=1, padx=5, pady=5)

        # Similar entries and dropdowns for other paths

        # Save button
        save_button = tk.Button(paths_window, text="Save and Close", command=paths_window.destroy)
        save_button.grid(row=1, column=0, columnspan=2, pady=10)

    def open_calibrate_window(self):
        # Example implementation for "Calibrate Defaults" window
        calibrate_window = tk.Toplevel(self.root)
        calibrate_window.title("Calibrate Defaults Configuration")

        # Add Entry widgets for each calibrate default
        calibrate_default1_label = tk.Label(calibrate_window, text="Calibrate Default 1:")
        calibrate_default1_label.grid(row=0, column=0, padx=5, pady=5)
        calibrate_default1_entry = tk.Entry(calibrate_window)
        calibrate_default1_entry.grid(row=0, column=1, padx=5, pady=5)

        # Similar entries for other calibrate defaults

        # Save button
        save_button = tk.Button(calibrate_window, text="Save and Close", command=calibrate_window.destroy)
        save_button.grid(row=1, column=0, columnspan=2, pady=10)

    def open_extract_window(self):
        # Example implementation for "Extract Defaults" window
        extract_window = tk.Toplevel(self.root)
        extract_window.title("Extract Defaults Configuration")

        # Add Entry widgets for each extract default
        extract_default1_label = tk.Label(extract_window, text="Extract Default 1:")
        extract_default1_label.grid(row=0, column=0, padx=5, pady=5)
        extract_default1_entry = tk.Entry(extract_window)
        extract_default1_entry.grid(row=0, column=1, padx=5, pady=5)

        # Similar entries for other extract defaults

        # Save button
        save_button = tk.Button(extract_window, text="Save and Close", command=extract_window.destroy)
        save_button.grid(row=1, column=0, columnspan=2, pady=10)

    def open_identify_window(self):
        # Example implementation for "Identify Defaults" window
        identify_window = tk.Toplevel(self.root)
        identify_window.title("Identify Defaults Configuration")

        # Add Entry widgets for each identify default
        identify_default1_label = tk.Label(identify_window, text="Identify Default 1:")
        identify_default1_label.grid(row=0, column=0, padx=5, pady=5)
        identify_default1_entry = tk.Entry(identify_window)
        identify_default1_entry.grid(row=0, column=1, padx=5, pady=5)

        # Similar entries for other identify defaults

        # Save button
        save_button = tk.Button(identify_window, text="Save and Close", command=identify_window.destroy)
        save_button.grid(row=1, column=0, columnspan=2, pady=10)

    def update_plots(self):
        # Example implementation for "Update Plots" button
        self.plot_rt()
        self.plot_mass()

    def plot_rt(self):
        # Example implementation for Retention Time Plot
        fig, ax = plt.subplots(figsize=(5, 3), tight_layout=True)
        ax.set_title("Retention Time Plot")
        ax.set_xlabel("Retention Time")
        ax.set_ylabel("Intensity")
        canvas = FigureCanvasTkAgg(fig, master=self.rt_plot_frame)
        canvas.draw()
        canvas.get_tk_widget().grid(row=1, column=0, sticky="nsew")
        self.rt_plot_frame.columnconfigure(0, weight=1)

    def plot_mass(self):
        # Example implementation for Mass Plot
        fig, ax = plt.subplots(figsize=(5, 3), tight_layout=True)
        ax.set_title("Mass Plot")
        ax.set_xlabel("Mass")
        ax.set_ylabel("Intensity")
        canvas = FigureCanvasTkAgg(fig, master=self.mass_plot_frame)
        canvas.draw()
        canvas.get_tk_widget().grid(row=1, column=0, sticky="nsew")
        self.mass_plot_frame.columnconfigure(0, weight=1)

if __name__ == "__main__":
    root = tk.Tk()
    app = MyWidgetApp(root)
    root.mainloop()






fig, ax = plt.subplots(figsize=(5, 3), tight_layout=True)


# before works

import tkinter as tk
from tkinter import ttk
from tkinter import messagebox
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import matplotlib.pyplot as plt

class MyWidgetApp:
    def __init__(self, root):
        self.root = root
        self.root.title("My Widget App")

        # Create main frame with three horizontal sections
        self.main_frame = ttk.Frame(root)
        self.main_frame.pack(expand=1, fill="both")

        # Upper left: Settings menu
        self.create_settings_menu()

        # Workflow section
        self.create_workflow_section()

        # Upper middle: Retention Time Plot
        self.rt_plot_frame = ttk.Frame(self.main_frame)
        self.rt_plot_frame.grid(row=1, column=1, sticky="nsew")
        self.plot_rt()

        # Lower middle: Mass Plot
        self.mass_plot_frame = ttk.Frame(self.main_frame)
        self.mass_plot_frame.grid(row=2, column=1, sticky="nsew")
        self.plot_mass()

        # Set row and column weights for resizing
        self.main_frame.rowconfigure(1, weight=1)
        self.main_frame.rowconfigure(2, weight=1)
        self.main_frame.columnconfigure(1, weight=1)

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
        defaults_menu.add_command(label="Identify", command=self.open_identify_window)

        # Close submenu
        settings_menu.add_command(label="Close", command=self.root.destroy)

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

        # Identification Section
        identification_section = ttk.Frame(self.workflow_frame, borderwidth=2, relief="groove")
        identification_section.grid(row=0, column=3, padx=5, pady=5)
        self.add_workflow_section(identification_section, "Identification", "do_alpinac", "do_matching")

    def add_workflow_section(self, frame, title, command1, command2=None):
        # Title label
        title_label = tk.Label(frame, text=title, font=("Helvetica", 10, "bold"))
        title_label.grid(row=0, column=0, pady=5)

        # Update Plot button
        if title == "Raw Data":
            update_plot_button = tk.Button(frame, text="Update Plots", command=self.update_plots)
            update_plot_button.grid(row=1, column=0, pady=5)
        elif title == "Mass Calibration":
            calibrate_button = tk.Button(frame, text="Calibrate", command=lambda: self.show_message(command1))
            calibrate_button.grid(row=1, column=0, pady=5)
        elif title == "Extraction":
            extract_button = tk.Button(frame, text="Extract", command=lambda: self.show_message(command1))
            extract_button.grid(row=1, column=0, pady=5)
        elif title == "Identification":
            alpinac_button = tk.Button(frame, text="Alpinac", command=lambda: self.show_message(command1))
            alpinac_button.grid(row=1, column=0, pady=5)
            matching_button = tk.Button(frame, text="Matching", command=lambda: self.show_message(command2))
            matching_button.grid(row=2, column=0, pady=5)

        # Message box
        self.create_message_box(frame)

    def create_message_box(self, frame):
        # Message box
        message_var = tk.StringVar()
        message_box = tk.Entry(frame, textvariable=message_var, state='readonly', width=20)
        message_box.grid(row=3, column=0, pady=5)

    def show_message(self, command):
        # Example implementation to update message box
        result = self.execute_workflow_command(command)
        message_var.set(result)

    def execute_workflow_command(self, command):
        # Example implementation of workflow commands
        if command == "do_import":
            return "Raw data imported successfully!"
        elif command == "do_calibration":
            return "Calibration successful!"
        elif command == "do_extraction":
            return "Extraction successful!"
        elif command == "do_alpinac":
            return "Alpinac successful!"
        elif command == "do_matching":
            return "Matching successful!"

    def open_paths_window(self):
        # Example implementation for "Paths" window
        paths_window = tk.Toplevel(self.root)
        paths_window.title("Paths Configuration")

        # Add Entry widgets and Dropdowns for each path
        path1_label = tk.Label(paths_window, text="Path 1:")
        path1_label.grid(row=0, column=0, padx=5, pady=5)
        path1_entry = tk.Entry(paths_window)
        path1_entry.grid(row=0, column=1, padx=5, pady=5)

        # Similar entries and dropdowns for other paths

        # Save button
        save_button = tk.Button(paths_window, text="Save and Close", command=paths_window.destroy)
        save_button.grid(row=1, column=0, columnspan=2, pady=10)

    def open_calibrate_window(self):
        # Example implementation for "Calibrate Defaults" window
        calibrate_window = tk.Toplevel(self.root)
        calibrate_window.title("Calibrate Defaults Configuration")

        # Add Entry widgets for each calibrate default
        calibrate_default1_label = tk.Label(calibrate_window, text="Calibrate Default 1:")
        calibrate_default1_label.grid(row=0, column=0, padx=5, pady=5)
        calibrate_default1_entry = tk.Entry(calibrate_window)
        calibrate_default1_entry.grid(row=0, column=1, padx=5, pady=5)

        # Similar entries for other calibrate defaults

        # Save button
        save_button = tk.Button(calibrate_window, text="Save and Close", command=calibrate_window.destroy)
        save_button.grid(row=1, column=0, columnspan=2, pady=10)

    def open_extract_window(self):
        # Example implementation for "Extract Defaults" window
        extract_window = tk.Toplevel(self.root)
        extract_window.title("Extract Defaults Configuration")

        # Add Entry widgets for each extract default
        extract_default1_label = tk.Label(extract_window, text="Extract Default 1:")
        extract_default1_label.grid(row=0, column=0, padx=5, pady=5)
        extract_default1_entry = tk.Entry(extract_window)
        extract_default1_entry.grid(row=0, column=1, padx=5, pady=5)

        # Similar entries for other extract defaults

        # Save button
        save_button = tk.Button(extract_window, text="Save and Close", command=extract_window.destroy)
        save_button.grid(row=1, column=0, columnspan=2, pady=10)

    def open_identify_window(self):
        # Example implementation for "Identify Defaults" window
        identify_window = tk.Toplevel(self.root)
        identify_window.title("Identify Defaults Configuration")

        # Add Entry widgets for each identify default
        identify_default1_label = tk.Label(identify_window, text="Identify Default 1:")
        identify_default1_label.grid(row=0, column=0, padx=5, pady=5)
        identify_default1_entry = tk.Entry(identify_window)
        identify_default1_entry.grid(row=0, column=1, padx=5, pady=5)

        # Similar entries for other identify defaults

        # Save button
        save_button = tk.Button(identify_window, text="Save and Close", command=identify_window.destroy)
        save_button.grid(row=1, column=0, columnspan=2, pady=10)

    def update_plots(self):
        # Example implementation for "Update Plots" button
        self.plot_rt()
        self.plot_mass()

    def plot_rt(self):
        # Example implementation for Retention Time Plot
        fig, ax = plt.subplots(figsize=(5, 3), tight_layout=True)
        ax.set_title("Retention Time Plot")
        ax.set_xlabel("Retention Time")
        ax.set_ylabel("Intensity")

        # Navigation toolbar
        toolbar = NavigationToolbar2Tk(canvas, self.rt_plot_frame)
        toolbar.update()
        canvas.get_tk_widget().grid(row=1, column=0, columnspan=2, sticky="nsew")

        # Entry fields for mass start and mass end
        mass_start_label = tk.Label(self.rt_plot_frame, text="Mass Start:")
        mass_start_label.grid(row=2, column=0, padx=5, pady=5)
        mass_start_entry = tk.Entry(self.rt_plot_frame)
        mass_start_entry.grid(row=2, column=1, padx=5, pady=5)

        mass_end_label = tk.Label(self.rt_plot_frame, text="Mass End:")
        mass_end_label.grid(row=3, column=0, padx=5, pady=5)
        mass_end_entry = tk.Entry(self.rt_plot_frame)
        mass_end_entry.grid(row=3, column=1, padx=5, pady=5)

        # Save button for mass range
        mass_save_button = tk.Button(self.rt_plot_frame, text="Apply Mass Range", command=lambda: self.apply_mass_range(ax, mass_start_entry.get(), mass_end_entry.get()))
        mass_save_button.grid(row=4, column=0, columnspan=2, pady=10)

        # Canvas for plotting
        canvas = FigureCanvasTkAgg(fig, master=self.rt_plot_frame)
        canvas.draw()
        canvas.get_tk_widget().grid(row=1, column=0, columnspan=2, sticky="nsew")
        self.rt_plot_frame.columnconfigure(0, weight=1)

    def apply_mass_range(self, ax, mass_start, mass_end):
        # Example implementation to apply mass range in the plot
        print("Applying Mass Range:", mass_start, mass_end)
        # Your logic to update the plot with the provided mass range

    def plot_mass(self):
        # Example implementation for Mass Plot
        fig, ax = plt.subplots(figsize=(5, 3), tight_layout=True)
        ax.set_title("Mass Plot")
        ax.set_xlabel("Mass")
        ax.set_ylabel("Intensity")

        # Navigation toolbar
        toolbar = NavigationToolbar2Tk(canvas, self.mass_plot_frame)
        toolbar.update()
        canvas.get_tk_widget().grid(row=1, column=0, columnspan=2, sticky="nsew")

        # Entry fields for RT start and RT end
        rt_start_label = tk.Label(self.mass_plot_frame, text="RT Start:")
        rt_start_label.grid(row=2, column=0, padx=5, pady=5)
        rt_start_entry = tk.Entry(self.mass_plot_frame)
        rt_start_entry.grid(row=2, column=1, padx=5, pady=5)

        rt_end_label = tk.Label(self.mass_plot_frame, text="RT End:")
        rt_end_label.grid(row=3, column=0, padx=5, pady=5)
        rt_end_entry = tk.Entry(self.mass_plot_frame)
        rt_end_entry.grid(row=3, column=1, padx=5, pady=5)

        # Save button for RT range
        rt_save_button = tk.Button(self.mass_plot_frame, text="Apply RT Range", command=lambda: self.apply_rt_range(ax, rt_start_entry.get(), rt_end_entry.get()))
        rt_save_button.grid(row=4, column=0, columnspan=2, pady=10)

        # Canvas for plotting
        canvas = FigureCanvasTkAgg(fig, master=self.mass_plot_frame)
        canvas.draw()
        canvas.get_tk_widget().grid(row=1, column=0, columnspan=2, sticky="nsew")
        self.mass_plot_frame.columnconfigure(0, weight=1)

    def apply_rt_range(self, ax, rt_start, rt_end):
        # Example implementation to apply RT range in the plot
        print("Applying RT Range:", rt_start, rt_end)
        # Your logic to update the plot with the provided RT range

if __name__ == "__main__":
    root = tk.Tk()
    app = MyWidgetApp(root)
    root.mainloop()




import tkinter as tk
from tkinter import ttk
from tkinter import messagebox
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import matplotlib.pyplot as plt

class MyWidgetApp:
    def __init__(self, root):
        self.root = root
        self.root.title("My Widget App")

        # Create main frame with three horizontal sections
        self.main_frame = ttk.Frame(root)
        self.main_frame.pack(expand=1, fill="both")

        # Upper left: Settings menu
        self.create_settings_menu()

        # Workflow section
        self.create_workflow_section()

        # Upper middle: Retention Time Plot
        self.rt_plot_frame = ttk.Frame(self.main_frame)
        self.rt_plot_frame.grid(row=1, column=1, sticky="nsew")
        self.plot_rt()

        # Lower middle: Mass Plot
        self.mass_plot_frame = ttk.Frame(self.main_frame)
        self.mass_plot_frame.grid(row=2, column=1, sticky="nsew")
        self.plot_mass()

        # Set row and column weights for resizing
        self.main_frame.rowconfigure(1, weight=1)
        self.main_frame.rowconfigure(2, weight=1)
        self.main_frame.columnconfigure(1, weight=1)

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
        defaults_menu.add_command(label="Identify", command=self.open_identify_window)

        # Close submenu
        settings_menu.add_command(label="Close", command=self.root.destroy)

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

        # Identification Section
        identification_section = ttk.Frame(self.workflow_frame, borderwidth=2, relief="groove")
        identification_section.grid(row=0, column=3, padx=5, pady=5)
        self.add_workflow_section(identification_section, "Identification", "do_alpinac", "do_matching")

    def add_workflow_section(self, frame, title, command1, command2=None):
        # Title label
        title_label = tk.Label(frame, text=title, font=("Helvetica", 10, "bold"))
        title_label.grid(row=0, column=0, pady=5)

        # Update Plot button
        if title == "Raw Data":
            update_plot_button = tk.Button(frame, text="Update Plots", command=self.update_plots)
            update_plot_button.grid(row=1, column=0, pady=5)
        elif title == "Mass Calibration":
            calibrate_button = tk.Button(frame, text="Calibrate", command=lambda: self.show_message(command1))
            calibrate_button.grid(row=1, column=0, pady=5)
        elif title == "Extraction":
            extract_button = tk.Button(frame, text="Extract", command=lambda: self.show_message(command1))
            extract_button.grid(row=1, column=0, pady=5)
        elif title == "Identification":
            alpinac_button = tk.Button(frame, text="Alpinac", command=lambda: self.show_message(command1))
            alpinac_button.grid(row=1, column=0, pady=5)
            matching_button = tk.Button(frame, text="Matching", command=lambda: self.show_message(command2))
            matching_button.grid(row=2, column=0, pady=5)

        # Message box
        self.create_message_box(frame)

    def create_message_box(self, frame):
        # Message box
        message_var = tk.StringVar()
        message_box = tk.Entry(frame, textvariable=message_var, state='readonly', width=20)
        message_box.grid(row=3, column=0, pady=5)

    def show_message(self, command):
        # Example implementation to update message box
        result = self.execute_workflow_command(command)
        message_var.set(result)

    def execute_workflow_command(self, command):
        # Example implementation of workflow commands
        if command == "do_import":
            return "Raw data imported successfully!"
        elif command == "do_calibration":
            return "Calibration successful!"
        elif command == "do_extraction":
            return "Extraction successful!"
        elif command == "do_alpinac":
            return "Alpinac successful!"
        elif command == "do_matching":
            return "Matching successful!"

    def open_paths_window(self):
        # Example implementation for "Paths" window
        paths_window = tk.Toplevel(self.root)
        paths_window.title("Paths Configuration")

        # Add Entry widgets and Dropdowns for each path
        path1_label = tk.Label(paths_window, text="Path 1:")
        path1_label.grid(row=0, column=0, padx=5, pady=5)
        path1_entry = tk.Entry(paths_window)
        path1_entry.grid(row=0, column=1, padx=5, pady=5)

        # Similar entries and dropdowns for other paths

        # Save button
        save_button = tk.Button(paths_window, text="Save and Close", command=paths_window.destroy)
        save_button.grid(row=1, column=0, columnspan=2, pady=10)

    def open_calibrate_window(self):
        # Example implementation for "Calibrate Defaults" window
        calibrate_window = tk.Toplevel(self.root)
        calibrate_window.title("Calibrate Defaults Configuration")

        # Add Entry widgets for each calibrate default
        calibrate_default1_label = tk.Label(calibrate_window, text="Calibrate Default 1:")
        calibrate_default1_label.grid(row=0, column=0, padx=5, pady=5)
        calibrate_default1_entry = tk.Entry(calibrate_window)
        calibrate_default1_entry.grid(row=0, column=1, padx=5, pady=5)

        # Similar entries for other calibrate defaults

        # Save button
        save_button = tk.Button(calibrate_window, text="Save and Close", command=calibrate_window.destroy)
        save_button.grid(row=1, column=0, columnspan=2, pady=10)

    def open_extract_window(self):
        # Example implementation for "Extract Defaults" window
        extract_window = tk.Toplevel(self.root)
        extract_window.title("Extract Defaults Configuration")

        # Add Entry widgets for each extract default
        extract_default1_label = tk.Label(extract_window, text="Extract Default 1:")
        extract_default1_label.grid(row=0, column=0, padx=5, pady=5)
        extract_default1_entry = tk.Entry(extract_window)
        extract_default1_entry.grid(row=0, column=1, padx=5, pady=5)

        # Similar entries for other extract defaults

        # Save button
        save_button = tk.Button(extract_window, text="Save and Close", command=extract_window.destroy)
        save_button.grid(row=1, column=0, columnspan=2, pady=10)

    def open_identify_window(self):
        # Example implementation for "Identify Defaults" window
        identify_window = tk.Toplevel(self.root)
        identify_window.title("Identify Defaults Configuration")

        # Add Entry widgets for each identify default
        identify_default1_label = tk.Label(identify_window, text="Identify Default 1:")
        identify_default1_label.grid(row=0, column=0, padx=5, pady=5)
        identify_default1_entry = tk.Entry(identify_window)
        identify_default1_entry.grid(row=0, column=1, padx=5, pady=5)

        # Similar entries for other identify defaults

        # Save button
        save_button = tk.Button(identify_window, text="Save and Close", command=identify_window.destroy)
        save_button.grid(row=1, column=0, columnspan=2, pady=10)

    def update_plots(self):
        # Example implementation for "Update Plots" button
        self.plot_rt()
        self.plot_mass()

    def plot_rt(self):
        # Example implementation for Retention Time Plot
        fig, ax = plt.subplots(figsize=(5, 3), tight_layout=True)
        ax.set_title("Retention Time Plot")
        ax.set_xlabel("Retention Time")
        ax.set_ylabel("Intensity")

        # Navigation toolbar
        toolbar = NavigationToolbar2Tk(canvas, self.rt_plot_frame)
        toolbar.update()
        canvas = FigureCanvasTkAgg(fig, master=self.rt_plot_frame)
        canvas.draw()
        canvas.get_tk_widget().grid(row=1, column=0, columnspan=2, sticky="nsew")
        self.rt_plot_frame.columnconfigure(0, weight=1)

    def plot_mass(self):
        # Example implementation for Mass Plot
        fig, ax = plt.subplots(figsize=(5, 3), tight_layout=True)
        ax.set_title("Mass Plot")
        ax.set_xlabel("Mass")
        ax.set_ylabel("Intensity")

        # Navigation toolbar
        toolbar = NavigationToolbar2Tk(canvas, self.mass_plot_frame)
        toolbar.update()
        canvas = FigureCanvasTkAgg(fig, master=self.mass_plot_frame)
        canvas.draw()
        canvas.get_tk_widget().grid(row=1, column=0, columnspan=2, sticky="nsew")
        self.mass_plot_frame.columnconfigure(0, weight=1)

    def apply_mass_range(self, ax, mass_start, mass_end):
        # Example implementation to apply mass range in the plot
        print("Applying Mass Range:", mass_start, mass_end)
        # Your logic to update the plot with the provided mass range

    def apply_rt_range(self, ax, rt_start, rt_end):
        # Example implementation to apply RT range in the plot
        print("Applying RT Range:", rt_start, rt_end)
        # Your logic to update the plot with the provided RT range

if __name__ == "__main__":
    root = tk.Tk()
    app = MyWidgetApp(root)
    root.mainloop()

# two previous not working

import tkinter as tk
from tkinter import ttk
from tkinter import messagebox
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import matplotlib.pyplot as plt

class MyWidgetApp:
    def __init__(self, root):
        self.root = root
        self.root.title("My Widget App")

        # Create main frame with three horizontal sections
        self.main_frame = ttk.Frame(root)
        self.main_frame.pack(expand=1, fill="both")

        # Upper left: Settings menu
        self.create_settings_menu()

        # Workflow section
        self.create_workflow_section()

        # Upper middle: Retention Time Plot
        self.rt_plot_frame = ttk.Frame(self.main_frame)
        self.rt_plot_frame.grid(row=1, column=1, sticky="nsew")
        self.plot_rt()

        # Lower middle: Mass Plot
        self.mass_plot_frame = ttk.Frame(self.main_frame)
        self.mass_plot_frame.grid(row=2, column=1, sticky="nsew")
        self.plot_mass()

        # Set row and column weights for resizing
        self.main_frame.rowconfigure(1, weight=1)
        self.main_frame.rowconfigure(2, weight=1)
        self.main_frame.columnconfigure(1, weight=1)

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
        defaults_menu.add_command(label="Identify", command=self.open_identify_window)

        # Close submenu
        settings_menu.add_command(label="Close", command=self.root.destroy)

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

        # Identification Section
        identification_section = ttk.Frame(self.workflow_frame, borderwidth=2, relief="groove")
        identification_section.grid(row=0, column=3, padx=5, pady=5)
        self.add_workflow_section(identification_section, "Identification", "do_alpinac", "do_matching")

    def add_workflow_section(self, frame, title, command1, command2=None):
        # Title label
        title_label = tk.Label(frame, text=title, font=("Helvetica", 10, "bold"))
        title_label.grid(row=0, column=0, pady=5)

        # Update Plot button
        if title == "Raw Data":
            update_plot_button = tk.Button(frame, text="Update Plots", command=self.update_plots)
            update_plot_button.grid(row=1, column=0, pady=5)
        elif title == "Mass Calibration":
            calibrate_button = tk.Button(frame, text="Calibrate", command=lambda: self.show_message(command1))
            calibrate_button.grid(row=1, column=0, pady=5)
        elif title == "Extraction":
            extract_button = tk.Button(frame, text="Extract", command=lambda: self.show_message(command1))
            extract_button.grid(row=1, column=0, pady=5)
        elif title == "Identification":
            alpinac_button = tk.Button(frame, text="Alpinac", command=lambda: self.show_message(command1))
            alpinac_button.grid(row=1, column=0, pady=5)
            matching_button = tk.Button(frame, text="Matching", command=lambda: self.show_message(command2))
            matching_button.grid(row=2, column=0, pady=5)

        # Message box
        self.create_message_box(frame)

    def create_message_box(self, frame):
        # Message box
        message_var = tk.StringVar()
        message_box = tk.Entry(frame, textvariable=message_var, state='readonly', width=20)
        message_box.grid(row=3, column=0, pady=5)

    def show_message(self, command):
        # Example implementation to update message box
        result = self.execute_workflow_command(command)
        message_var.set(result)

    def execute_workflow_command(self, command):
        # Example implementation of workflow commands
        if command == "do_import":
            return "Raw data imported successfully!"
        elif command == "do_calibration":
            return "Calibration successful!"
        elif command == "do_extraction":
            return "Extraction successful!"
        elif command == "do_alpinac":
            return "Alpinac successful!"
        elif command == "do_matching":
            return "Matching successful!"

    def open_paths_window(self):
        # Example implementation for "Paths" window
        paths_window = tk.Toplevel(self.root)
        paths_window.title("Paths Configuration")

        # Add Entry widgets and Dropdowns for each path
        path1_label = tk.Label(paths_window, text="Path 1:")
        path1_label.grid(row=0, column=0, padx=5, pady=5)
        path1_entry = tk.Entry(paths_window)
        path1_entry.grid(row=0, column=1, padx=5, pady=5)

        # Similar entries and dropdowns for other paths

        # Save button
        save_button = tk.Button(paths_window, text="Save and Close", command=paths_window.destroy)
        save_button.grid(row=1, column=0, columnspan=2, pady=10)

    def open_calibrate_window(self):
        # Example implementation for "Calibrate Defaults" window
        calibrate_window = tk.Toplevel(self.root)
        calibrate_window.title("Calibrate Defaults Configuration")

        # Add Entry widgets for each calibrate default
        calibrate_default1_label = tk.Label(calibrate_window, text="Calibrate Default 1:")
        calibrate_default1_label.grid(row=0, column=0, padx=5, pady=5)
        calibrate_default1_entry = tk.Entry(calibrate_window)
        calibrate_default1_entry.grid(row=0, column=1, padx=5, pady=5)

        # Similar entries for other calibrate defaults

        # Save button
        save_button = tk.Button(calibrate_window, text="Save and Close", command=calibrate_window.destroy)
        save_button.grid(row=1, column=0, columnspan=2, pady=10)

    def open_extract_window(self):
        # Example implementation for "Extract Defaults" window
        extract_window = tk.Toplevel(self.root)
        extract_window.title("Extract Defaults Configuration")

        # Add Entry widgets for each extract default
        extract_default1_label = tk.Label(extract_window, text="Extract Default 1:")
        extract_default1_label.grid(row=0, column=0, padx=5, pady=5)
        extract_default1_entry = tk.Entry(extract_window)
        extract_default1_entry.grid(row=0, column=1, padx=5, pady=5)

        # Similar entries for other extract defaults

        # Save button
        save_button = tk.Button(extract_window, text="Save and Close", command=extract_window.destroy)
        save_button.grid(row=1, column=0, columnspan=2, pady=10)

    def open_identify_window(self):
        # Example implementation for "Identify Defaults" window
        identify_window = tk.Toplevel(self.root)
        identify_window.title("Identify Defaults Configuration")

        # Add Entry widgets for each identify default
        identify_default1_label = tk.Label(identify_window, text="Identify Default 1:")
        identify_default1_label.grid(row=0, column=0, padx=5, pady=5)
        identify_default1_entry = tk.Entry(identify_window)
        identify_default1_entry.grid(row=0, column=1, padx=5, pady=5)

        # Similar entries for other identify defaults

        # Save button
        save_button = tk.Button(identify_window, text="Save and Close", command=identify_window.destroy)
        save_button.grid(row=1, column=0, columnspan=2, pady=10)

    def update_plots(self):
        # Example implementation for "Update Plots" button
        self.plot_rt()
        self.plot_mass()

    def plot_rt(self):
        # Example implementation for Retention Time Plot
        fig, ax = plt.subplots(figsize=(5, 3), tight_layout=True)
        ax.set_title("Retention Time Plot")
        ax.set_xlabel("Retention Time")
        ax.set_ylabel("Intensity")

        # Navigation toolbar
        toolbar = NavigationToolbar2Tk(canvas, self.rt_plot_frame)
        toolbar.update()
        canvas = FigureCanvasTkAgg(fig, master=self.rt_plot_frame)
        canvas.draw()
        canvas.get_tk_widget().grid(row=1, column=0, columnspan=2, sticky="nsew")
        self.rt_plot_frame.columnconfigure(0, weight=1)

    def plot_mass(self):
        # Example implementation for Mass Plot
        fig, ax = plt.subplots(figsize=(5, 3), tight_layout=True)
        ax.set_title("Mass Plot")
        ax.set_xlabel("Mass")
        ax.set_ylabel("Intensity")

        # Navigation toolbar
        toolbar = NavigationToolbar2Tk(canvas, self.mass_plot_frame)
        toolbar.update()
        canvas = FigureCanvasTkAgg(fig, master=self.mass_plot_frame)
        canvas.draw()
        canvas.get_tk_widget().grid(row=1, column=0, columnspan=2, sticky="nsew")
        self.mass_plot_frame.columnconfigure(0, weight=1)

    def apply_mass_range(self, ax, mass_start, mass_end):
        # Example implementation to apply mass range in the plot
        print("Applying Mass Range:", mass_start, mass_end)
        # Your logic to update the plot with the provided mass range

    def apply_rt_range(self, ax, rt_start, rt_end):
        # Example implementation to apply RT range in the plot
        print("Applying RT Range:", rt_start, rt_end)
        # Your logic to update the plot with the provided RT range

if __name__ == "__main__":
    root = tk.Tk()
    app = MyWidgetApp(root)
    root.mainloop()
