# The chain (Rumors reference intended)

# ----------------------------------------- Higher-level Chain functions -----------------------------------------------


def makeChain(functionDictList, preiterable=None, inMethod=None, in_kwargs=None, writeInfoExists=False,
              recordMethod=None, record_kwargs=None, writeInfo_recordkw=None, extractMethod=None, extract_kwargs=None):
    """
    FUNCTION: produces an object (class: Chain) that contains information about a source of data (preiterable) and, if
        required, a method by which the data can be converted to a iterable form (inMethod, in_kwargs). The object
        contains info for applying a series of functions, with specified parameters (functionDictList) at each
        iteration. Additionally, the method records information with each iteration which may be a list of objects
        returned from the final function in the data processing chain or some information pertaining to this depending
        on whether a function to be applied to this output is provided (recordMethod, record_kwargs,
        writeInfo_recordkw). This added functionality is useful in the case where one functionDictList is useful in
        the context of different types of data, requiring different types of records [apologies, this is vague but
        meaningful]. Once iterations are complete, the record will either be returned as is or operated on depending on
        whether a function has been supplied (extractMethod, extract_kwargs).
    ARGUMENTS:
        preiterable: data source (type depends on inMethod, if inMethod is None, this must be iterable)
        functionDictList: list of dictionaries, which must contain specific keys including: function, suffix,
            writeOut, args*, kwargs* (* optional). The value at function should be the function to execute.
            The value at writeOut should be None or a dictionary with keys function, args*, kwargs*, priorLink_kw*,
            newLink_kw, and suffix_kw*. ([{}, ...])
        inMethod: If the data source must be converted to an iterable object, a function can be supplied
        in_kwargs: If the inMethod function requires additional arguments, supply these here. (dict)
        writeInfoExists: Should data be written out for each iteration? (bool)
        recordMethod: function which can be applied to the final output of processing chain each iteration
            in order to produce an object which can be added to the record  (function or None).
            If None, the the final object of the processing chain will be added to the record.
        record_kwargs: key word arguments for the recordMethod function. (dict)
        writeInfo_recordkw: name of the argument to which any writeInfo (as above) should be supplied in the
            recordMethod function.
        extractMethod: function that can be applied to the accumulated record following the final iteration
            to produce the final result returned by Chain().execute() (function or None).
            If None, Chain().execute() will return the record as is.
        extract_kwargs: key word arguments for the extractMethod function (dict).
    DEPENDENCIES: Class: Chain
    RETURNS: instance of Chain
    """
    chain = Chain(preiterable, inMethod=inMethod, in_kwargs=in_kwargs, writeInfoExists=writeInfoExists,
                  recordMethod=recordMethod, record_kwargs=record_kwargs, writeInfo_recordkw=writeInfo_recordkw,
                  extractMethod=extractMethod, extract_kwargs=extract_kwargs, functionDictList=functionDictList)
    return chain


def convertChain(chain, preiterable, inMethod=None, in_kwargs=None, writeInfo=True,
                 recordMethod=None, extractMethod=None, writeInfo_recordkw=None):
    """
    FUNCTION: Alter an existing chain object in order to apply its chain functions to a different type of data
    """
    chain.newInput(preiterable)
    chain.newInputMethods(inMethod=inMethod, in_kwargs=in_kwargs, writeInfo=writeInfo)
    chain.newOutputMethods(recordMethod=recordMethod, extractMethod=extractMethod,
                           writeInfo_recordkw=writeInfo_recordkw)
    return chain


def alterChainMethod(chain, suffix, update, paramName=None):
    if paramName is None:
        chain.updateChainMethod(suffix, update)
    if paramName is not None:
        chain.updateChainMethodParam(suffix, paramName, update)
    return chain


# ---------------------------------------- Chain-related helper functions ----------------------------------------------

def makeChainTemplateDict(functionDictList=None, noLinks=0, writeOutBools=None):
    if functionDictList is None:
        functionDictList = []
        for link in range(noLinks):
            functionDict = makeFunctionDictTemplate(writeOut=writeOutBools[link])
            functionDictList.append(functionDict)

    chainDict = {'preiterable': None, 'in_kwargs': None, 'writeInfoExists': False,
                 'functionDictList': functionDictList,
                 'recordMethod': None, 'record_kwargs': None, 'writeInfo_recordkw': None,
                 'extractMethod': None, 'extract_kwargs': None}
    return chainDict


def makeFunctionDictTemplate(writeOut=False):
    writeOutDict = None
    if writeOut:
        writeOutDict = {'function': None, 'args': None, 'kwargs': None, 'priorLink_kw': None, 'newLink_kw': None,
                        'writeInfo_kw': None, 'suffix_kw': None}
    functionDict = {'function': None, 'suffix': None, 'args': None, 'kwargs': None, 'writeOut': writeOutDict}
    return functionDict


# ------------------------------------------------ Chain classes -------------------------------------------------------

class Chain:
    def __init__(self, preiterable=None,
                 inMethod=None,
                 in_kwargs=None,
                 writeInfoExists=False,
                 functionDictList=None,
                 recordMethod=None,
                 record_kwargs=None,
                 writeInfo_recordkw=None,
                 extractMethod=None,
                 extract_kwargs=None,
                 **kwargs):
        """
        DATA ITERATION
            :param preiterable: data source (type depends on inMethod and functions in chain). Not necessary and wont
                be used if data is supplied to the execute method.
            :param inMethod: method for converting data source into iterable object. If None, the preiterable object
                must be iterable (function)
            :param in_kwargs: key word arguments for inMethod function (dict)
            :param writeInfoExists: When iterating over the  inputted data (once inMethod(preiterable, **in_kwargs)
                is applied), does the data have the form ~ (chainObject, writeInfo)*?
                    * where chainObject is the data object with which to apply chain functions and
                    writeInfo is some information required to write out information at this step.
        CHAIN FUNCTIONS
            :param functionDictList:  list of dictionaries, which must contain specific keys including: function,
                suffix, writeOut, args*, kwargs* (* optional). The value at function should be the function to execute.
                The value at writeOut should be None or a dictionary with keys function, args*, kwargs*, priorLink_kw*,
                newLink_kw, and suffix_kw*. ([{}, ...])
        DATA RECORDING OVER ITERATIONS
            :param recordMethod: function which can be applied to the final output of processing chain each iteration
                in order to produce an object which can be added to the record (function or None).
                If None, the the final object of the processing chain will be added to the record.
            :param record_kwargs: key word arguments for the recordMethod function (dict).
            :param writeInfo_recordkw: name of the argument to which any writeInfo (as above) should be supplied in
                the recordMethod function.
            :param extractMethod: function that can be applied to the accumulated record following the final iteration
                to produce the final result returned by Chain().execute() (function or None).
                If None, Chain().execute() will return the record as is.
            :param extract_kwargs: key word arguments for the extractMethod function (dict).
        """
        # the following relate to the method by which data source is converted to an iterable object
        # and the
        self.preiterable = preiterable  # probably a substack  or directory of files
        self.inMethod = inMethod  # function for reading input into chain (e.g., generateFromDirectory(...))
        self.in_kwargs = in_kwargs
        if self.in_kwargs is None:
            self.in_kwargs = {}
        self.writeInfoExists = writeInfoExists
        # The following code defines the chain if functionDictList is specified
        self.chain = []  # list of ChainMethods
        if functionDictList is not None:
            self.defineChain(functionDictList=functionDictList)
        # define the scribe (can be overwritten)
        self.scribe = Scribe(recordMethod=recordMethod,
                             record_kwargs=record_kwargs,
                             writeInfo_recordkw=writeInfo_recordkw,
                             extractMethod=extractMethod,
                             extract_kwargs=extract_kwargs)

    # Basic methods ---------------------------------------
    def defineChain(self, functionDictList):
        for functionDict in functionDictList:
            newLink = ChainMethod(functionDict)
            self.chain.append(newLink)

    def addScribe(self, recordMethod=None, record_kwargs=None, writeInfo_recordkw=None,
                  extractMethod=None, extract_kwargs=None):
        self.scribe = Scribe(recordMethod=recordMethod,
                             record_kwargs=record_kwargs,
                             writeInfo_recordkw=writeInfo_recordkw,
                             extractMethod=extractMethod,
                             extract_kwargs=extract_kwargs)

    def addChainMethod(self, functionDict):
        newLink = ChainMethod(functionDict)
        self.chain.append(newLink)

    def removeChainMethod(self, suffix):
        for i in range(len(self.chain)):
            if self.chain[i].suffix == suffix:
                remove = i
        if remove is not None:
            del self.chain[remove]
        else:
            print('No ChainMethod with suffix {}'.format(suffix))

    def getChainMethodSuffixes(self):
        suffixes = [chainMethod.suffix for chainMethod in chain]
        return suffixes

    # Update methods ---------------------------------------
    def updateChainMethod(self, suffix, functionDict):
        """
        Redefine one instance of ChainMethod with the desired value at ChainMethod.suffix
        """
        for i in range(len(self.chain)):
            if self.chain[i].suffix == suffix:
                self.chain[i].update(functionDict)

    def updateChainMethodParam(self, suffix, paramName, update):
        """
        Redefine one parameter in an existing ChainMethod with the desired value at ChainMethod.suffix
        """
        for i in range(len(self.chain)):
            if self.chain[i].suffix == suffix:
                self.chain[i].funct[paramName] = update
                self.chain[i].setAll()

    def newInput(self, newInput):
        self.preiterable = newInput

    def newInputMethods(self, inMethod=None, in_kwargs=None, write=True):
        self.inMethod = inMethod
        self.in_kwargs = in_kwargs
        self.writeInfoExists = writeInfoExists

    def newOutputMethods(self, recordMethod=None, record_kwargs=None, extractMethod=None, extract_kwargs=None,
                         writeInfo_recordkw=None):
        self.scribe = Scribe(recordMethod=recordMethod,
                             record_kwargs=record_kwargs,
                             writeInfo_recordkw=writeInfo_recordkw,
                             extractMethod=extractMethod,
                             extract_kwargs=extract_kwargs)

    # Usage method -----------------------------------------
    def execute(self, data=None):
        # define the objects that must be iterated over
        if data is None:
            data = self.preiterable
        if self.writeInfoExists and self.inMethod is not None:
            for chainObj, writeInfo in self.inMethod(data, **self.in_kwargs):
                for chainMethod in self.chain:
                    chainObj = chainMethod.run(chainObj, writeInfo=writeInfo)
                self.scribe.addRecord(chainObj, writeInfo=writeInfo)
            output = self.scribe.extract()
            self.scribe.cleanRecord()
            return output
        if not self.writeInfoExists and self.inMethod is not None:
            for chainObj in self.inMethod(data, **self.in_kwargs):
                for chainMethod in self.chain:
                    chainObj = chainMethod.run(chainObj)
                self.scribe.addRecord(chainObj)
            output = self.scribe.extract()
            self.scribe.cleanRecord()
            return output
        else:
            for chainObj in data:
                for chainMethod in self.chain:
                    chainObj = chainMethod.run(chainObj)
                self.scribe.addRecord(chainObj)
            output = self.scribe.extract()
            self.scribe.cleanRecord()
            return output


class ChainMethod:
    """
    Takes a function dictionary as the initial argument. An instance of ChainMethod is an object that allows
        a function to be executed with a series of pre-defined parameters (as args or kwargs) when supplied only
        a single argument.
    If a writeOut key is supplied with an appropriate dict, data will be saved at this step.
        If a writeOut function is supplied, several methods for write out can be applied at this step. If additional
        information required for write out to occur, the run method should be supplied with a writeInfo argument of
        not None (likely str). This will be assumed to be the third argument of the function if present.
        If the dict contains a usePrior key of True, the main function's input argument (chainObj) will be used
        as the first argument of the writeOut function and the main function's output will be used as the second.
    """

    def __init__(self, functionDict):
        self.funct = functionDict
        self.setAll()

    # Basic methods ----------------------------------------
    def setAll(self):
        self.args = self.funct.get('args')
        if self.args is None:
            self.args = []
        self.kwargs = self.funct.get('kwargs')
        if self.kwargs is None:
            self.kwargs = {}
        self.writeOut = self.funct.get('writeOut')
        self.suffix = self.funct.get('suffix')
        self.writeInfo_kw = self.funct.get('writeInfo_kw')

    def function(self, chainObj, *args, **kwargs):
        funct = self.funct['function']
        return funct(chainObj, *args, **kwargs)

    # Update methods ---------------------------------------
    def update_args(self, update):
        self.args = update

    def update_kwargs(self, update):
        self.kwargs = update

    def update_suffix(self, update):
        self.suffix = update

    def update_writeInfo_kw(self, update):
        self.writeInfo_kw = update

    def update(self, update):
        self.funct = update
        self.setAll()

    # Usage methods ----------------------------------------
    def run(self, chainObj, writeInfo=None):
        """
        :param chainObj: the current data object on which to operate
        :param writeInfo: some information associated with said object
        :return: object to be passed to the next chain method or to the scribe
        """
        kwargs = self.kwargs.copy()
        if self.writeInfo_kw is not None:
            kwargs[self.writeInfo_kw] = writeInfo
        chainOut = self.function(chainObj, *self.args, **kwargs)
        if self.writeOut is not None:
            self.writeOutData(chainObj, chainOut, writeInfo=writeInfo)
        return chainOut

    def writeOutData(self, chainObj, chainOut, writeInfo=None):
        """
        :param chainObj: former data object (priorLink)
        :param chainOut: current data object (newLink) produced by self.function
        :param writeInfo: information required to write the file (str or None)
        :return: NA. The write function is expected to save one or more files.
        """
        writeFunc = self.writeOut['function']
        argsList = self.writeOut.get('args')
        if argsList is None:  # cannot be none internally. TypeError
            argsList = []
        kwargsDict = self.writeOut.get('kwargs').copy()
        if kwargsDict is None:  # cannot be none internally. TypeError
            kwargsDict = {}
        chainObj_kw = self.writeOut.get('priorLink_kw')
        chainOut_kw = self.writeOut.get('newLink_kw')
        writeInfo_kw = self.writeOut.get('writeInfo_kw')
        suffix_kw = self.writeOut.get('suffix_kw')
        kwargsDict[chainObj_kw] = chainObj
        kwargsDict[chainOut_kw] = chainOut
        if writeInfo is not None:
            kwargsDict[writeInfo_kw] = writeInfo
        if suffix_kw is not None:
            kwargsDict[self.suffix] = suffix_kw
        writeFunc(*argsList, **kwargsDict)


class Scribe:
    """
    This class contains an attribute record, which is a list of objects written out over an iteration. The method for
    recording the objects is inputted as recordMethod, which may also be None (in this case, the input of add record
    is simply added to the record list). The extractMethod if not None is a function applied to the record.
    """

    def __init__(self, recordMethod=None, record_kwargs=None, writeInfo_recordkw=None,
                 extractMethod=None, extract_kwargs=None):
        """
        :param recordMethod: function for producing objects to be added to the record attribute (which is a list)
        :param record_kwargs: key word arguments to be supplied to the recordMethod
        :param writeInfo_recordkw: key word (str) for argument of recordMethod that should be supplied with the
            writeInfo (e.g., a file name when processing a directory of files).
        :param extractMethod: function
        :param extract_kwargs:
        """
        self.record = []
        self.recordMethod = recordMethod
        self.record_kwargs = record_kwargs
        if self.record_kwargs is None:
            self.record_kwargs = {}  # prevents TypeError when self.addRecord is called + avoids mutable default 
        self.extractMethod = extractMethod
        self.extract_kwargs = extract_kwargs
        if self.extract_kwargs is None:
            self.extract_kwargs = {}  # prevents TypeError when self.extract is called
        self.writeInfo_recordkw = writeInfo_recordkw

    def addRecord(self, chainObj, writeInfo=None):
        if self.recordMethod is not None:
            if writeInfo is not None:
                self.record_kwargs[self.writeInfo_recordkw] = writeInfo
            rec = self.recordMethod(chainObj, **self.record_kwargs)
        else:
            rec = chainObj
        if rec is not None:
            self.record.append(rec)

    def extract(self):
        if self.extractMethod is not None:
            output = self.extractMethod(self.record, **self.extract_kwargs)
        else:
            output = self.record
        return output

    def cleanRecord(self):
        self.record = []
