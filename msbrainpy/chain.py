# The chain (Rumors reference intended)

# ----------------------------------------- Higher-level Chain functions -----------------------------------------------


def make_chain(function_dict_list, preiterable=None, in_method=None, in_kwargs=None, write_info_exists=False,
              record_method=None, record_kwargs=None, write_info_recordkw=None, extract_method=None, extract_kwargs=None):
    """
    FUNCTION: produces an object (class: Chain) that contains information about a source of data (preiterable) and, if
        required, a method by which the data can be converted to a iterable form (in_method, in_kwargs). The object
        contains info for applying a series of functions, with specified parameters (function_dict_list) at each
        iteration. Additionally, the method records information with each iteration which may be a list of objects
        returned from the final function in the data processing chain or some information pertaining to this depending
        on whether a function to be applied to this output is provided (record_method, record_kwargs,
        write_info_recordkw). This added functionality is useful in the case where one function_dict_list is useful in
        the context of different types of data, requiring different types of records [apologies, this is vague but
        meaningful]. Once iterations are complete, the record will either be returned as is or operated on depending on
        whether a function has been supplied (extract_method, extract_kwargs).
    ARGUMENTS:
        preiterable: data source (type depends on in_method, if in_method is None, this must be iterable)
        function_dict_list: list of dictionaries, which must contain specific keys including: function, suffix,
            write_out, args*, kwargs* (* optional). The value at function should be the function to execute.
            The value at write_out should be None or a dictionary with keys function, args*, kwargs*, prior_link_kw*,
            new_link_kw, and suffix_kw*. ([{}, ...])
        in_method: If the data source must be converted to an iterable object, a function can be supplied
        in_kwargs: If the in_method function requires additional arguments, supply these here. (dict)
        write_info_exists: Should data be written out for each iteration? (bool)
        record_method: function which can be applied to the final output of processing chain each iteration
            in order to produce an object which can be added to the record  (function or None).
            If None, the the final object of the processing chain will be added to the record.
        record_kwargs: key word arguments for the record_method function. (dict)
        write_info_recordkw: name of the argument to which any write_info (as above) should be supplied in the
            record_method function.
        extract_method: function that can be applied to the accumulated record following the final iteration
            to produce the final result returned by Chain().execute() (function or None).
            If None, Chain().execute() will return the record as is.
        extract_kwargs: key word arguments for the extract_method function (dict).
    DEPENDENCIES: Class: Chain
    RETURNS: instance of Chain
    """
    chain = Chain(preiterable, in_method=in_method, in_kwargs=in_kwargs, write_info_exists=write_info_exists,
                  record_method=record_method, record_kwargs=record_kwargs, write_info_recordkw=write_info_recordkw,
                  extract_method=extract_method, extract_kwargs=extract_kwargs, function_dict_list=function_dict_list)
    return chain


def convert_chain(chain, preiterable, in_method=None, in_kwargs=None, write_info=True,
                 record_method=None, extract_method=None, write_info_recordkw=None):
    """
    FUNCTION: Alter an existing chain object in order to apply its chain functions to a different type of data
    """
    chain.new_input(preiterable)
    chain.new_input_methods(in_method=in_method, in_kwargs=in_kwargs, write_info=write_info)
    chain.new_output_methods(record_method=record_method, extract_method=extract_method,
                           write_info_recordkw=write_info_recordkw)
    return chain


def alter_chain_method(chain, suffix, update, param_name=None):
    if param_name is None:
        chain.update_chain_method(suffix, update)
    if param_name is not None:
        chain.update_chain_method_param(suffix, param_name, update)
    return chain


# ---------------------------------------- Chain-related helper functions ----------------------------------------------

def make_chain_template_dict(function_dict_list=None, no_links=0, write_out_bools=None):
    if function_dict_list is None:
        function_dict_list = []
        for link in range(no_links):
            function_dict = make_function_dict_template(write_out=write_outBools[link])
            function_dict_list.append(function_dict)

    chain_dict = {'preiterable': None, 'in_kwargs': None, 'write_info_exists': False,
                 'function_dict_list': function_dict_list,
                 'record_method': None, 'record_kwargs': None, 'write_info_recordkw': None,
                 'extract_method': None, 'extract_kwargs': None}
    return chain_dict


def make_function_dict_template(write_out=False):
    write_out_dict = None
    if write_out:
        write_out_dict = {'function': None, 'args': None, 'kwargs': None, 'prior_link_kw': None, 'new_link_kw': None,
                        'write_info_kw': None, 'suffix_kw': None}
    function_dict = {'function': None, 'suffix': None, 'args': None, 'kwargs': None, 'write_out': write_outDict}
    return function_dict


# ------------------------------------------------ Chain classes -------------------------------------------------------

class Chain:
    def __init__(self, preiterable=None,
                 in_method=None,
                 in_kwargs=None,
                 write_info_exists=False,
                 function_dict_list=None,
                 record_method=None,
                 record_kwargs=None,
                 write_info_recordkw=None,
                 extract_method=None,
                 extract_kwargs=None,
                 **kwargs):
        """
        DATA ITERATION
            :param preiterable: data source (type depends on in_method and functions in chain). Not necessary and wont
                be used if data is supplied to the execute method.
            :param in_method: method for converting data source into iterable object. If None, the preiterable object
                must be iterable (function)
            :param in_kwargs: key word arguments for in_method function (dict)
            :param write_info_exists: When iterating over the  inputted data (once in_method(preiterable, **in_kwargs)
                is applied), does the data have the form ~ (chain_object, write_info)*?
                    * where chain_object is the data object with which to apply chain functions and
                    write_info is some information required to write out information at this step.
        CHAIN FUNCTIONS
            :param function_dict_list:  list of dictionaries, which must contain specific keys including: function,
                suffix, write_out, args*, kwargs* (* optional). The value at function should be the function to execute.
                The value at write_out should be None or a dictionary with keys function, args*, kwargs*, prior_link_kw*,
                new_link_kw, and suffix_kw*. ([{}, ...])
        DATA RECORDING OVER ITERATIONS
            :param record_method: function which can be applied to the final output of processing chain each iteration
                in order to produce an object which can be added to the record (function or None).
                If None, the the final object of the processing chain will be added to the record.
            :param record_kwargs: key word arguments for the record_method function (dict).
            :param write_info_recordkw: name of the argument to which any write_info (as above) should be supplied in
                the record_method function.
            :param extract_method: function that can be applied to the accumulated record following the final iteration
                to produce the final result returned by Chain().execute() (function or None).
                If None, Chain().execute() will return the record as is.
            :param extract_kwargs: key word arguments for the extract_method function (dict).
        """
        # the following relate to the method by which data source is converted to an iterable object
        # and the
        self.preiterable = preiterable  # probably a substack  or directory of files
        self.in_method = in_method  # function for reading input into chain (e.g., generate_from_directory(...))
        self.in_kwargs = in_kwargs
        if self.in_kwargs is None:
            self.in_kwargs = {}
        self.write_info_exists = write_info_exists
        # The following code defines the chain if function_dict_list is specified
        self.chain = []  # list of chain_methods
        if function_dict_list is not None:
            self.define_chain(function_dict_list=function_dict_list)
        # define the scribe (can be overwritten)
        self.scribe = Scribe(record_method=record_method,
                             record_kwargs=record_kwargs,
                             write_info_recordkw=write_info_recordkw,
                             extract_method=extract_method,
                             extract_kwargs=extract_kwargs)

    # Basic methods ---------------------------------------
    def define_chain(self, function_dict_list):
        for function_dict in function_dictList:
            new_link = chain_method(function_dict)
            self.chain.append(new_link)

    def add_scribe(self, record_method=None, record_kwargs=None, write_info_recordkw=None,
                  extract_method=None, extract_kwargs=None):
        self.scribe = Scribe(record_method=record_method,
                             record_kwargs=record_kwargs,
                             write_info_recordkw=write_info_recordkw,
                             extract_method=extract_method,
                             extract_kwargs=extract_kwargs)

    def add_chain_method(self, function_dict):
        new_link = chain_method(function_dict)
        self.chain.append(new_link)

    def remove_chain_method(self, suffix):
        for i in range(len(self.chain)):
            if self.chain[i].suffix == suffix:
                remove = i
        if remove is not None:
            del self.chain[remove]
        else:
            print('No chain_method with suffix {}'.format(suffix))

    def get_chain_method_suffixes(self):
        suffixes = [chain_method.suffix for chain_method in chain]
        return suffixes

    # Update methods ---------------------------------------
    def update_chain_method(self, suffix, function_dict):
        """
        Redefine one instance of chain_method with the desired value at chain_method.suffix
        """
        for i in range(len(self.chain)):
            if self.chain[i].suffix == suffix:
                self.chain[i].update(function_dict)

    def update_chain_method_param(self, suffix, param_name, update):
        """
        Redefine one parameter in an existing chain_method with the desired value at chain_method.suffix
        """
        for i in range(len(self.chain)):
            if self.chain[i].suffix == suffix:
                self.chain[i].funct[param_name] = update
                self.chain[i].set_all()

    def new_input(self, new_input):
        self.preiterable = new_input

    def new_input_methods(self, in_method=None, in_kwargs=None, write=True):
        self.in_method = in_method
        self.in_kwargs = in_kwargs
        self.write_info_exists = write_info_exists

    def new_output_methods(self, record_method=None, record_kwargs=None, extract_method=None, extract_kwargs=None,
                         write_info_recordkw=None):
        self.scribe = Scribe(record_method=record_method,
                             record_kwargs=record_kwargs,
                             write_info_recordkw=write_info_recordkw,
                             extract_method=extract_method,
                             extract_kwargs=extract_kwargs)

    # Usage method -----------------------------------------
    def execute(self, data=None):
        # define the objects that must be iterated over
        if data is None:
            data = self.preiterable
        if self.write_info_exists and self.in_method is not None:
            for chain_obj, write_info in self.in_method(data, **self.in_kwargs):
                for chain_method in self.chain:
                    chain_obj = chain_method.run(chain_obj, write_info=write_info)
                self.scribe.add_record(chain_obj, write_info=write_info)
            output = self.scribe.extract()
            self.scribe.clean_record()
            return output
        if not self.write_info_exists and self.in_method is not None:
            for chain_obj in self.in_method(data, **self.in_kwargs):
                for chain_method in self.chain:
                    chain_obj = chain_method.run(chain_obj)
                self.scribe.add_record(chain_obj)
            output = self.scribe.extract()
            self.scribe.clean_record()
            return output
        else:
            for chain_obj in data:
                for chain_method in self.chain:
                    chain_obj = chain_method.run(chain_obj)
                self.scribe.add_record(chain_obj)
            output = self.scribe.extract()
            self.scribe.clean_record()
            return output


class ChainMethod:
    """
    Takes a function dictionary as the initial argument. An instance of chain_method is an object that allows
        a function to be executed with a series of pre-defined parameters (as args or kwargs) when supplied only
        a single argument.
    If a write_out key is supplied with an appropriate dict, data will be saved at this step.
        If a write_out function is supplied, several methods for write out can be applied at this step. If additional
        information required for write out to occur, the run method should be supplied with a write_info argument of
        not None (likely str). This will be assumed to be the third argument of the function if present.
        If the dict contains a use_prior key of True, the main function's input argument (chain_obj) will be used
        as the first argument of the write_out function and the main function's output will be used as the second.
    """

    def __init__(self, function_dict):
        self.funct = function_dict
        self.set_all()

    # Basic methods ----------------------------------------
    def set_all(self):
        self.args = self.funct.get('args')
        if self.args is None:
            self.args = []
        self.kwargs = self.funct.get('kwargs')
        if self.kwargs is None:
            self.kwargs = {}
        self.write_out = self.funct.get('write_out')
        self.suffix = self.funct.get('suffix')
        self.write_info_kw = self.funct.get('write_info_kw')

    def function(self, chain_obj, *args, **kwargs):
        funct = self.funct['function']
        return funct(chain_obj, *args, **kwargs)

    # Update methods ---------------------------------------
    def update_args(self, update):
        self.args = update

    def update_kwargs(self, update):
        self.kwargs = update

    def update_suffix(self, update):
        self.suffix = update

    def update_write_info_kw(self, update):
        self.write_info_kw = update

    def update(self, update):
        self.funct = update
        self.set_all()

    # Usage methods ----------------------------------------
    def run(self, chain_obj, write_info=None):
        """
        :param chain_obj: the current data object on which to operate
        :param write_info: some information associated with said object
        :return: object to be passed to the next chain method or to the scribe
        """
        kwargs = self.kwargs.copy()
        if self.write_info_kw is not None:
            kwargs[self.write_info_kw] = write_info
        chain_out = self.function(chain_obj, *self.args, **kwargs)
        if self.write_out is not None:
            self.write_out_data(chain_obj, chain_out, write_info=write_info)
        return chain_out

    def write_out_data(self, chain_obj, chain_out, write_info=None):
        """
        :param chain_obj: former data object (prior_link)
        :param chain_out: current data object (new_link) produced by self.function
        :param write_info: information required to write the file (str or None)
        :return: NA. The write function is expected to save one or more files.
        """
        write_func = self.write_out['function']
        args_list = self.write_out.get('args')
        if args_list is None:  # cannot be none internally. type_error
            args_list = []
        kwargs_dict = self.write_out.get('kwargs').copy()
        if kwargs_dict is None:  # cannot be none internally. type_error
            kwargs_dict = {}
        chain_obj_kw = self.write_out.get('prior_link_kw')
        chain_out_kw = self.write_out.get('new_link_kw')
        write_info_kw = self.write_out.get('write_info_kw')
        suffix_kw = self.write_out.get('suffix_kw')
        kwargs_dict[chain_obj_kw] = chain_obj
        kwargs_dict[chain_out_kw] = chain_out
        if write_info is not None:
            kwargs_dict[write_info_kw] = write_info
        if suffix_kw is not None:
            kwargs_dict[self.suffix] = suffix_kw
        write_func(*args_list, **kwargs_dict)


class Scribe:
    """
    This class contains an attribute record, which is a list of objects written out over an iteration. The method for
    recording the objects is inputted as record_method, which may also be None (in this case, the input of add record
    is simply added to the record list). The extract_method if not None is a function applied to the record.
    """

    def __init__(self, record_method=None, record_kwargs=None, write_info_recordkw=None,
                 extract_method=None, extract_kwargs=None):
        """
        :param record_method: function for producing objects to be added to the record attribute (which is a list)
        :param record_kwargs: key word arguments to be supplied to the record_method
        :param write_info_recordkw: key word (str) for argument of record_method that should be supplied with the
            write_info (e.g., a file name when processing a directory of files).
        :param extract_method: function
        :param extract_kwargs:
        """
        self.record = []
        self.record_method = record_method
        self.record_kwargs = record_kwargs
        if self.record_kwargs is None:
            self.record_kwargs = {}  # prevents type_error when self.add_record is called + avoids mutable default 
        self.extract_method = extract_method
        self.extract_kwargs = extract_kwargs
        if self.extract_kwargs is None:
            self.extract_kwargs = {}  # prevents type_error when self.extract is called
        self.write_info_recordkw = write_info_recordkw

    def add_record(self, chain_obj, write_info=None):
        if self.record_method is not None:
            if write_info is not None:
                self.record_kwargs[self.write_info_recordkw] = write_info
            rec = self.record_method(chain_obj, **self.record_kwargs)
        else:
            rec = chain_obj
        if rec is not None:
            self.record.append(rec)

    def extract(self):
        if self.extract_method is not None:
            output = self.extract_method(self.record, **self.extract_kwargs)
        else:
            output = self.record
        return output

    def clean_record(self):
        self.record = []
