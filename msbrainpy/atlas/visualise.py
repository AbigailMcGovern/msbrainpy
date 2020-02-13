import matplotlib.pyplot as plt
from msbrainpy.amba.wrangle import matchRows, amalgamateDFs

# -------------------------------------------------- Globals -----------------------------------------------------------
defaultColours = {'prefrontal': 'purple', 'somatomotor': 'pink', 'anterolateral': 'blue', 'temporal': 'green',
                  'visual': 'yellow', 'medial': 'red', 'undefined': 'black'}


# --------------------------------------------- Region scatter plots ---------------------------------------------------

def plotRegions(pointsArr, df, acronymHeader='acronym', groupHeader='group', name=None, xlab='x label',
                ylab='y label', size=(12, 12), dpi=300, colours=defaultColours):
    """
    FUNCTION: scater plot of brain regions with labled acronyms
    ARGUMENTS:
        pointsArr = np.ndarray, 2 variables with n rows of data, shape (n, 2)
        df = pd.DataFrame containing information about one of the variables (same len, correct order)
        acronymHeader = header of column in df containing region acronyms with which to annotate points (str)
        groupHeader = header of column in df containing information about groups by which to colour code points (str)
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
    acronyms = [acronym for acronym in df[acronymHeader]]
    fig, ax = plt.subplots()
    for i in range(len(pointsArr)):
        x, y = pointsArr[i]
        ax.annotate(acronyms[i], (x, y))
    ax.scatter(pointsArr[:, 0], pointsArr[:, 1], c=df[groupHeader].apply(lambda x: colours[x]))
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    fig.set_size_inches(size)
    plt.show()
    if name is not None:
        fig.savefig(name, dpi=dpi)
    return fig, ax


def plotRegions_Layers(pointsArr, df, name=None, acronymHeader='acronym', groupHeader='group', size=(14, 14),
                       dpi=300, colours=defaltColours):
    """
    FUNCTION: Plot values pertaining to layers of cortical regions with only layer annotated.
    ARGUMENTS:
        pointsArr = np.ndarray, 2 variables with n rows of data, shape (n, 2)
        df = pd.DataFrame containing information about one of the variables (same len, correct order)
        acronymHeader = header of column in df containing region acronyms with which to annotate points (str)
        groupHeader = header of column in df containing information about groups by which to colour code points (str)
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
    acronyms = [acronym for acronym in df[acronymHeader]]
    fig, ax = plt.subplots()
    for i in range(len(df)):
        x, y = pointsArr[i]
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
    ax.scatter(pointsArr[:, 0], pointsArr[:, 1], c=df[groupHeader].apply(lambda x: colours[x]))
    plt.show()
    if name is not None:
        fig.set_size_inches(size)
        fig.savefig(name, dpi=dpi)
    return fig, ax


def getPoints(df0, df1, values=('expression_density', 'expression_density'),
              groups=('structure_acronym', 'structure_acronym'), saveName=None, df_prefixes=None):
    """
    FUNCTION: obtains an array of xy points (2 variables with n rows of points, shape (n, 2)) from 2 data frames
        in order to obtain a labeled & colour-coded scatter plot (e.g., labels = acronyms, colours = groups).
    ARGUMENTS:
        df0_ = pd.DataFrame from which to obtain column 0 (x or independent)
        df1_ = pd.DataFrame from which to obtain column 1 (y or dependent)
        values = names of columns from which to obtain x and y, respectively (string, string)
        groups = names of columns containing the label assigned to these points (string, string)
        saveName = None OR name with which to save an amalgamated data frame
    DEPENDENCIES: Packages/modules/etc: numpy as np, pandas as pd
        Native: (1) msbrainpy.wrangle.amalgamateDFs(df_name_dict)
        1) makes a conjoint data frame from several smaller data frames with matching rows
    RETURNS: np.ndarray
    """
    df0, df1 = matchRows(df0, df1, groups)
    if saveName is not None:
        df = amalgamateDFs({df_prefixes[0]: df0, df_prefixes[1]: df1})
        df.to_csv(saveName)
    pointsArr = np.array([df0[values[0]], df1[values[1]]]).T
    return pointsArr


# ---------------------------------------- And now... plotting lines ---------------------------------------------------

def plotModel(df, colours, lines, namesHeader='names', valuesHeader='values',
              xlab='pvalb expression density', ylab='plp expression density',
              max_=0.16, min_=0.0, step=0.01):
    values = [i for i in df[valuesHeader]]
    names = [i for i in df[namesHeader]]
    fig, ax = plt.subplots()
    for i in range(len(values)):
        value = values[i]
        colour = colours[names[i]]
        line = lines[names[i]]
        x, y = getLinePoints(value[0], value[1], max_=max_, min_=max_, step=step)
        ax.plot(x, y, color=colour, linestyle=line)
        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)
    plt.legend(names)


def getLinePoints(intercept, slope, max_=0.16, min_=0.0, step=0.03):
    x = np.arange(0.0, 0.16, 0.01)  # np.arange(max_, min_, step)
    y = intercept + slope * x
    return x, y


# --------------------------------------------- Only vaguely related  --------------------------------------------------

def saveDictAsJSON(name, obj):
    """
    Save a dictionary containing a useful mapping of categories to display settings
    (e.g., groups to colours or line types)
    """
    json_ = json.dumps(obj)
    with open(name, "w") as file:
        file.write(json_)

