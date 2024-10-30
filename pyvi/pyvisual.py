import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

from typing import List, Dict

import re

class PyVi:
    class PyViSection:
        def __init__(self, name, param):
            self.name = name
            self.param = param
            self.iterations = {}

        def parse_line(self, iter, vector):
            self.iterations[int(iter)] = np.fromstring(vector, dtype=np.float64, sep=',')
        
        def plot_iter(self, param, k, color, label_axes, fmt):
            plt.plot(param, self.iterations[k], fmt, color=color)
            if label_axes:
                plt.xlabel(self.param)
                plt.ylabel(self.name)

    def __load_from_file(self, file):
        with open(file) as f:
            lines = f.readlines()

            GLOBAL_MODE = 0
            PARAMETER_MODE = 1
            SECTION_MODE = 2
            SECTION_MODE_VEC = 3

            state = GLOBAL_MODE

            current_section = ''

            for line in lines:
                if m := re.match(r'\[(.*?)\]', line):
                    if m.group(1) == 'Parameter':    
                        state = PARAMETER_MODE
                    elif m.group(1) == 'Section':
                        state = SECTION_MODE
                    continue
                
                if m := re.match(r'(.*?):(.*)', line):
                    self.params[m.group(1)] = np.fromstring(m.group(2), dtype=np.float64, sep=',')
                
                if m := re.match(r'\((.*?)\)\-\>\[(.*?)\]', line):
                    self.sections[m.group(1)] = PyVi.PyViSection(m.group(1), m.group(2))
                    current_section = m.group(1)
                    state = SECTION_MODE_VEC
                
                if m := re.match(r'I\[(\d+)\]\=(.*)', line):
                    self.sections[current_section].parse_line(m.group(1), m.group(2))

    def __init__(self, filename):
        self.params : Dict[str, np.ndarray] = {}
        self.sections : Dict[str, PyVi.PyViSection] = {}
        self.__load_from_file(filename)
    
    def plot_iteration(self, section_name, iteration_no, label_axes=False, plt_format='o', color=None):
        param = self.params[self.sections[section_name].param]
        self.sections[section_name].plot_iter(param, iteration_no, color, label_axes, plt_format)

    # plot all sections that depend on parameter
    def plot_common(self, parameter : str, iter_no : int, plt_format='o'):
        i = 0
        for key, section in self.sections.items():
            if section.param != parameter:
                continue
            ccolor = mcolors.TABLEAU_COLORS[list(mcolors.TABLEAU_COLORS.keys())[i]]
            plt.plot(self.params[parameter], section.iterations[iter_no], color=ccolor, label=key)
            i = (i + 1) % len(mcolors.TABLEAU_COLORS)
        plt.xlabel(parameter)
        plt.legend()

    # plot these sections, that depend on the same parameter
    def plot_sections(self, parameter : str, sections : List[str], iter_no):
        i = 0
        for key in sections:
            section = self.sections[key]
            if section.param != parameter:
                continue
            ccolor = mcolors.TABLEAU_COLORS[list(mcolors.TABLEAU_COLORS.keys())[i]]
            plt.plot(self.params[parameter], section.iterations[iter_no], color=ccolor, label=key)
            i = (i + 1) % len(mcolors.TABLEAU_COLORS)
        plt.xlabel(parameter)
        plt.legend()


