import numpy as np
import matplotlib.pyplot as plt
from msbrainpy.atlas.wrangle import match_rows, amalgamate_dFs

# -------------------------------------------------- Globals -----------------------------------------------------------
default_colours = {'prefrontal': 'purple', 'somatomotor': 'pink', 'anterolateral': 'blue', 'temporal': 'green',
                  'visual': 'yellow', 'medial': 'red', 'undefined': 'black'}


# --------------------------------------------- Region scatter plots ---------------------------------------------------

def plot_regions(points_arr, df, acronym_header='acronym', group_header='group', name=None, xlab='x label',
                ylab='y label', size=(12, 12), dpi=300, colours=default_colours):
    """
    FUNCTION: scater plot of brain regions with labled acronyms
    ARGUMENTS:
        points_arr = np.ndarray, 2 variables with n rows of data, shape (n, 2)
        df = pd.data_frame containing information about one of the variables (same len, correct order)
        acronym_header = header of column in df containing region acronyms with which to annotate points (str)
        group_header = header of column in df containing information about groups by which to colour code points (str)
        name = name with which to save file - i.e., some_name.png (str)
        xlab = name of x axis (str)
        ylab = name of y axis (str)
        size = size of figure in inches ((int, int))
        dpi = dots per inch (int)
        colours = dictionary. keys: categories from the group column, values: colours to assign points
            default -> {'prefrontal': 'purple', 'somatomotor': 'pink', 'anterolateral': 'blue', 'temporal': 'green',
                 'visual': 'yellow', 'medial': 'red', 'undefined': 'black'}
    DEPENDENCIES: Packages/modules/etc: matplotlib.pyplot as plt
    RETURNS: fig and ax objects
    Note: large image size required for densely populated scatter plots.
        12 x 12 inches may not be sufficient (just enough spread for all cortical layers)
    """
    acronyms = [acronym for acronym in df[acronym_header]]
    fig, ax = plt.subplots()
    for i in range(len(points_arr)):
        x, y = points_arr[i]
        ax.annotate(acronyms[i], (x, y))
    ax.scatter(points_arr[:, 0], points_arr[:, 1], c=df[group_header].apply(lambda x: colours[x]))
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    fig.set_size_inches(size)
    plt.show()
    if name is not None:
        fig.savefig(name, dpi=dpi)
    return fig, ax


def plot_regions__layers(points_arr, df, name=None, acronym_header='acronym', group_header='group', size=(14, 14),
                       dpi=300, colours=default_colours):
    """
    FUNCTION: Plot values pertaining to layers of cortical regions with only layer annotated.
    ARGUMENTS:
        points_arr = np.ndarray, 2 variables with n rows of data, shape (n, 2)
        df = pd.data_frame containing information about one of the variables (same len, correct order)
        acronym_header = header of column in df containing region acronyms with which to annotate points (str)
        group_header = header of column in df containing information about groups by which to colour code points (str)
        name = name with which to save file - i.e., some_name.png (str)
        xlab = name of x axis (str)
        ylab = name of y axis (str)
        size = size of figure in inches ((int, int))
        dpi = dots per inch (int)
        colours = dictionary. keys: categories from the group column, values: colours to assign points
    DEPENDENCIES: Packages/modules/etc: matplotlib.pyplot as plt
    RETURNS: fig and ax objects

    Note: large image size required for densely populated scatter plots.
        12 x 12 inches may not be sufficient (just enough spread for all cortical layers)
    """
    acronyms = [acronym for acronym in df[acronym_header]]
    fig, ax = plt.subplots()
    for i in range(len(df)):
        x, y = points_arr[i]
        if '2' in acronyms[i]:
            if 'FRP' in acronyms[i]:
                ax.annotate(acronyms[i], (x, y))
            if 'FRP' not in acronyms[i]:
                ax.annotate('2/3', (x, y))
        if '1' in acronyms[i]:
            ax.annotate('1', (x, y))
        if '4' in acronyms[i]:
            ax.annotate('4', (x, y))
        if '5' in acronyms[i]:
            ax.annotate('5', (x, y))
        if '6' in acronyms[i]:
            ax.annotate('6', (x, y))
    ax.scatter(points_arr[:, 0], points_arr[:, 1], c=df[group_header].apply(lambda x: colours[x]))
    plt.show()
    if name is not None:
        fig.set_size_inches(size)
        fig.savefig(name, dpi=dpi)
    return fig, ax


def get_points(df0, df1, values=('expression_density', 'expression_density'),
              groups=('structure_acronym', 'structure_acronym'), save_name=None, df_prefixes=None):
    """
    FUNCTION: obtains an array of xy points (2 variables with n rows of points, shape (n, 2)) from 2 data frames
        in order to obtain a labeled & colour-coded scatter plot (e.g., labels = acronyms, colours = groups).
    ARGUMENTS:
        df0_ = pd.data_frame from which to obtain column 0 (x or independent)
        df1_ = pd.data_frame from which to obtain column 1 (y or dependent)
        values = names of columns from which to obtain x and y, respectively (string, string)
        groups = names of columns containing the label assigned to these points (string, string)
        save_name = None OR name with which to save an amalgamated data frame
    DEPENDENCIES: Packages/modules/etc: numpy as np, pandas as pd
        Native: (1) msbrainpy.wrangle.amalgamate_d_fs(df_name_dict)
        1) makes a conjoint data frame from several smaller data frames with matching rows
    RETURNS: np.ndarray
    """
    df0, df1 = match_rows(df0, df1, groups)
    if save_name is not None:
        df = amalgamate_d_fs({df_prefixes[0]: df0, df_prefixes[1]: df1})
        df.to_csv(save_name)
    points_arr = np.array([df0[values[0]], df1[values[1]]]).T
    return points_arr


# ---------------------------------------- And now... plotting lines ---------------------------------------------------

def plot_model(df, colours, lines, names_header='names', values_header='values',
              xlab='pvalb expression density', ylab='plp expression density',
              max_=0.16, min_=0.0, step=0.01):
    values = [i for i in df[values_header]]
    names = [i for i in df[names_header]]
    fig, ax = plt.subplots()
    for i in range(len(values)):
        value = values[i]
        colour = colours[names[i]]
        line = lines[names[i]]
        x, y = get_line_points(value[0], value[1], max_=max_, min_=max_, step=step)
        ax.plot(x, y, color=colour, linestyle=line)
        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)
    plt.legend(names)


def get_line_points(intercept, slope, max_=0.16, min_=0.0, step=0.03):
    x = np.arange(0.0, 0.16, 0.01)  # np.arange(max_, min_, step)
    y = intercept + slope * x
    return x, y


# --------------------------------------------- Only vaguely related  --------------------------------------------------

def save_dict_as_json(name, obj):
    """
    Save a dictionary containing a useful mapping of categories to display settings
    (e.g., groups to colours or line types)
    """
    json_ = json.dumps(obj)
    with open(name, "w") as file:
        file.write(json_)

