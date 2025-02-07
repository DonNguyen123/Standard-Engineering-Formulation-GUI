import tkinter as tk
from tkinter import ttk
import importlib
import inspect
import os
import sys
from pathlib import Path
import ast
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

class LibraryInterface:
    def __init__(self, root):
        self.root = root
        self.root.title("SEF Library Interface")
        self.root.geometry("1000x600")  # Made wider to accommodate the documentation panel
        
        # Create main container with two columns
        self.main_frame = ttk.Frame(root, padding="10")
        self.main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # Configure grid weights
        self.main_frame.columnconfigure(1, weight=1)  # Make documentation column expand
        
        # Left panel for existing controls
        self.left_panel = ttk.Frame(self.main_frame)
        self.left_panel.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # Right panel for documentation
        self.create_documentation_panel()
        
        # Create frames in left panel
        self.create_library_frame()
        self.create_function_frame()
        self.create_parameter_frame()
        self.create_output_frame()
        
        # Load libraries
        self.libraries = {}
        self.load_libraries()
        
    def create_library_frame(self):
        # Library selection frame
        library_frame = ttk.LabelFrame(self.left_panel, text="Select Library", padding="5")
        library_frame.grid(row=0, column=0, padx=5, pady=(0,2), sticky=(tk.W, tk.E))
        
        self.library_var = tk.StringVar()
        self.library_dropdown = ttk.Combobox(library_frame, textvariable=self.library_var, width=40)  # Increased width
        self.library_dropdown.grid(row=0, column=0, sticky=(tk.W, tk.E))
        self.library_dropdown.bind('<<ComboboxSelected>>', self.on_library_select)

    def create_function_frame(self):
        # Function selection frame
        function_frame = ttk.LabelFrame(self.left_panel, text="Select Function", padding="5")
        function_frame.grid(row=1, column=0, padx=5, pady=2, sticky=(tk.W, tk.E))
        
        self.function_var = tk.StringVar()
        self.function_dropdown = ttk.Combobox(function_frame, textvariable=self.function_var, width=40)  # Increased width
        self.function_dropdown.grid(row=0, column=0, sticky=(tk.W, tk.E))
        self.function_dropdown.bind('<<ComboboxSelected>>', self.on_function_select)
        
    def create_parameter_frame(self):
        # Parameter input frame
        self.parameter_frame = ttk.LabelFrame(self.left_panel, text="Input Parameters", padding="5")
        self.parameter_frame.grid(row=2, column=0, padx=5, pady=2, sticky=(tk.W, tk.E))
        
        # Will be populated dynamically based on function selection
        self.parameter_entries = {}

    def create_output_frame(self):
        # Output frame
        output_frame = ttk.LabelFrame(self.left_panel, text="Output", padding="5")
        output_frame.grid(row=3, column=0, padx=5, pady=2, sticky=(tk.W, tk.E))
        
        # Create a notebook for tabs
        self.output_notebook = ttk.Notebook(output_frame)
        self.output_notebook.grid(row=0, column=0, padx=5, pady=5, sticky=(tk.W, tk.E))
        
        # Text output tab
        self.text_frame = ttk.Frame(self.output_notebook)
        self.output_notebook.add(self.text_frame, text='Data')
        
        self.output_text = tk.Text(self.text_frame, height=20, width=60)
        self.output_text.grid(row=0, column=0, padx=5, pady=5)
        
        # Figure output tab
        self.figure_frame = ttk.Frame(self.output_notebook)
        self.output_notebook.add(self.figure_frame, text='Figure')
        
        # Zoom control frame
        self.zoom_frame = ttk.Frame(self.figure_frame)
        self.zoom_frame.grid(row=0, column=0, sticky='ew', padx=5, pady=2)
        
        # Zoom dropdown
        self.zoom_var = tk.StringVar(value='100%')
        self.zoom_dropdown = ttk.Combobox(self.zoom_frame, 
                                        textvariable=self.zoom_var,
                                        values=['50%', '100%', '200%'],
                                        width=10,
                                        state='readonly')
        self.zoom_dropdown.grid(row=0, column=0, padx=5)
        self.zoom_dropdown.bind('<<ComboboxSelected>>', self.on_zoom_change)
        
        # Create a canvas with scrollbar for the figure
        self.canvas_container = ttk.Frame(self.figure_frame)
        self.canvas_container.grid(row=1, column=0, sticky='nsew')
        
        # Make the figure frame expandable
        self.figure_frame.grid_rowconfigure(1, weight=1)
        self.figure_frame.grid_columnconfigure(0, weight=1)
        
        # Create canvas and scrollbar
        self.fig_canvas = tk.Canvas(self.canvas_container)
        self.scrollbar = ttk.Scrollbar(self.canvas_container, orient="vertical", command=self.fig_canvas.yview)
        
        # Configure canvas scrolling
        self.fig_canvas.configure(yscrollcommand=self.scrollbar.set)
        
        # Grid canvas and scrollbar
        self.fig_canvas.grid(row=0, column=0, sticky='nsew')
        self.scrollbar.grid(row=0, column=1, sticky='ns')
        
        # Make the canvas container expandable
        self.canvas_container.grid_rowconfigure(0, weight=1)
        self.canvas_container.grid_columnconfigure(0, weight=1)
        
        # Calculate button
        self.calculate_button = ttk.Button(output_frame, text="Calculate", command=self.calculate)
        self.calculate_button.grid(row=1, column=0, pady=5)

    def on_zoom_change(self, event=None):
        if hasattr(self, 'current_figure'):
            self.display_figure(self.current_figure)

    def display_figure(self, fig):
        # Store current figure for zoom operations
        self.current_figure = fig
        
        # Clear previous figure
        for widget in self.fig_canvas.winfo_children():
            widget.destroy()
        
        # Get zoom level
        zoom = int(self.zoom_var.get().rstrip('%'))
        scale_factor = zoom / 100.0
        
        # Get original figure size
        fig_width, fig_height = fig.get_size_inches()
        dpi = fig.get_dpi()
        
        # Calculate new size
        new_width = int(fig_width * dpi * scale_factor)
        new_height = int(fig_height * dpi * scale_factor)
        
        # Configure canvas size
        self.fig_canvas.configure(width=new_width)
        
        # Display the figure
        figure_canvas = FigureCanvasTkAgg(fig, self.fig_canvas)
        figure_canvas.draw()
        
        # Get the figure widget and add it to the canvas
        figure_widget = figure_canvas.get_tk_widget()
        figure_widget.configure(width=new_width, height=new_height)
        
        # Add figure to canvas
        window = self.fig_canvas.create_window((0, 0), window=figure_widget, anchor='nw')
        
        # Update scroll region after the figure is added
        self.fig_canvas.update_idletasks()
        self.fig_canvas.configure(scrollregion=self.fig_canvas.bbox("all"))
        
    def create_documentation_panel(self):
        # Documentation frame
        doc_frame = ttk.LabelFrame(self.main_frame, text="Function Documentation", padding="5")
        doc_frame.grid(row=0, column=1, padx=5, pady=5, sticky=(tk.N, tk.S, tk.E, tk.W))
        
        self.doc_text = tk.Text(doc_frame, width=50, wrap=tk.WORD)
        self.doc_text.grid(row=0, column=0, sticky=(tk.N, tk.S, tk.E, tk.W))
        
        # Add scrollbar
        doc_scroll = ttk.Scrollbar(doc_frame, orient="vertical", command=self.doc_text.yview)
        doc_scroll.grid(row=0, column=1, sticky=(tk.N, tk.S))
        self.doc_text.configure(yscrollcommand=doc_scroll.set)
        
    def load_libraries(self):
        
        def get_application_path():
            if getattr(sys, 'frozen', False):
                # If the application is run as a bundle (exe)
                application_path = sys._MEIPASS
            else:
                # If the application is run from a Python interpreter
                application_path = os.path.dirname(os.path.abspath(__file__))
            return application_path
        
        # Get the application path
        current_dir = Path(get_application_path())
        current_file = Path(__file__).name.split('.')[0] if not getattr(sys, 'frozen', False) else "your_main_script_name"
        
        # Add the directory to Python's path
        if current_dir not in sys.path:
            sys.path.insert(0, str(current_dir))
            
        try:
            python_files = [f.stem for f in current_dir.glob('*.py') 
                        if f.suffix == '.py' 
                        and not f.stem.startswith('__')
                        and f.stem != current_file]
            
            for module_name in python_files:
                try:
                    # Create full path to module
                    module_path = current_dir / f"{module_name}.py"
                    
                    # Create spec and register module
                    spec = importlib.util.spec_from_file_location(module_name, module_path)
                    if spec:
                        module = importlib.util.module_from_spec(spec)
                        # Register the module in sys.modules before executing
                        sys.modules[module_name] = module
                        spec.loader.exec_module(module)
                        self.libraries[module_name] = module
                    else:
                        print(f"Could not create spec for {module_name}")
                except Exception as e:
                    print(f"Error importing {module_name}: {e}")
                    # Clean up sys.modules if import failed
                    if module_name in sys.modules:
                        del sys.modules[module_name]
                    
            # Update library dropdown
            self.library_dropdown['values'] = list(self.libraries.keys())
            
        except Exception as e:
            print(f"Error loading libraries: {e}")
            
    def on_library_select(self, event):
        selected_library = self.library_var.get()
        if selected_library in self.libraries:
            module = self.libraries[selected_library]
            
            # Get only the functions that are defined in this module
            # Filter out imported functions by checking their module origin
            functions = [name for name, obj in inspect.getmembers(module)
                        if inspect.isfunction(obj) and 
                        obj.__module__ == selected_library and  # Only get functions from this module
                        not name.startswith('_')]
            
            # Update function dropdown
            self.function_dropdown['values'] = functions
            self.function_var.set('')  # Clear current selection
            
            # Clear parameter frame
            for widget in self.parameter_frame.winfo_children():
                widget.destroy()
            self.parameter_entries.clear()
            
    def on_function_select(self, event):
        selected_library = self.library_var.get()
        selected_function = self.function_var.get()
        
        if selected_library in self.libraries and selected_function:
            module = self.libraries[selected_library]
            function = getattr(module, selected_function)
            
            # Clear existing parameter widgets
            for widget in self.parameter_frame.winfo_children():
                widget.destroy()
            self.parameter_entries.clear()
            
            # Get function parameters
            params = inspect.signature(function).parameters
            
            # Create input fields for each parameter
            for i, (param_name, param) in enumerate(params.items()):
                label = ttk.Label(self.parameter_frame, text=f"{param_name}:")
                label.grid(row=i, column=0, padx=5, pady=2, sticky=tk.W)
                
                entry = ttk.Entry(self.parameter_frame)
                entry.grid(row=i, column=1, padx=5, pady=2, sticky=(tk.W, tk.E))
                
                # If parameter has a default value, show it as placeholder
                if param.default != inspect.Parameter.empty:
                    entry.insert(0, str(param.default))
                
                self.parameter_entries[param_name] = entry
            
            # Update documentation
            self.update_documentation(function)
            
    def get_function_docstring(self, module, function_name):
        """
        Extract docstring that appears above the function in source code,
        handling different file encodings.
        """
        try:
            # Get the source file path
            source_file = inspect.getsourcefile(module)
            
            # Try different encodings
            encodings = ['utf-8', 'iso-8859-1', 'cp1252', 'latin1']
            module_source = None
            
            for encoding in encodings:
                try:
                    with open(source_file, 'r', encoding=encoding) as f:
                        module_source = f.read()
                    break  # If successful, break the loop
                except UnicodeDecodeError:
                    continue
                    
            if module_source is None:
                print("Could not read source file with any supported encoding")
                return None
                
            # Parse the source code
            tree = ast.parse(module_source)
            
            # Get line numbers for all string literals and functions
            doc_candidates = {}  # Store all string literals and their line numbers
            target_function_line = None
            
            for node in ast.walk(tree):
                if isinstance(node, ast.FunctionDef) and node.name == function_name:
                    target_function_line = node.lineno
                elif isinstance(node, ast.Expr) and isinstance(node.value, ast.Str):
                    doc_candidates[node.lineno] = node.value.s
            
            if target_function_line is None:
                return None
                
            # Find the closest docstring that appears before the function
            closest_doc_line = None
            closest_docstring = None
            
            for doc_line, docstring in doc_candidates.items():
                if doc_line < target_function_line:
                    if closest_doc_line is None or doc_line > closest_doc_line:
                        closest_doc_line = doc_line
                        closest_docstring = docstring
            
            return closest_docstring
                        
        except Exception as e:
            print(f"Error extracting docstring: {e}")
            return None

    def update_documentation(self, function):
        # Clear current documentation
        self.doc_text.delete(1.0, tk.END)
        
        # Get the module containing the function
        module = sys.modules[function.__module__]
        
        # Try to get the docstring that appears above the function
        doc = self.get_function_docstring(module, function.__name__)
        
        if doc:
            self.doc_text.insert(tk.END, doc)
        else:
            # Fallback to regular docstring if no above-docstring found
            doc = function.__doc__
            if doc:
                self.doc_text.insert(tk.END, doc)
            else:
                self.doc_text.insert(tk.END, "No documentation available for this function.")
                
    def calculate(self):
        selected_library = self.library_var.get()
        selected_function = self.function_var.get()
        
        if selected_library in self.libraries and selected_function:
            module = self.libraries[selected_library]
            function = getattr(module, selected_function)
            
            try:
                params = {}
                for param_name, entry in self.parameter_entries.items():
                    value = entry.get()
                    try:
                        value = ast.literal_eval(value)
                    except (ValueError, SyntaxError):
                        pass
                    params[param_name] = value
                
                # Call function and get results
                result = function(**params)
                
                # Clear previous outputs
                self.output_text.delete(1.0, tk.END)
                
                if isinstance(result, tuple) and len(result) == 2:
                    fig, data = result
                    
                    # Display the figure with current zoom level
                    self.display_figure(fig)
                    
                    # Display the data
                    self.output_text.insert(tk.END, "Results:\n\n")
                    for key, value in data.items():
                        self.output_text.insert(tk.END, f"{key}:\n{value}\n\n")
                else:
                    self.output_text.insert(tk.END, f"Result: {result}")
                
            except Exception as e:
                self.output_text.delete(1.0, tk.END)
                self.output_text.insert(tk.END, f"Error: {str(e)}")

def main():
    root = tk.Tk()
    app = LibraryInterface(root)
    root.mainloop()

if __name__ == "__main__":
    main()