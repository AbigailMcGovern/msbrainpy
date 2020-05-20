import os
import re


# rename functions and variables from camelCase to snake_case (weep for the cleanliness of it)
def get_relevant_directories(directory, exclude_patterns=['.git', 'scripts']):
    walk_generator = os.walk(directory)
    directory_list = []
    for iterable in walk_generator:
        path = iterable[0]
        for exclude_pattern in exclude_patterns:
            if path.find(exclude_pattern) == -1:
                if path.find('setup.py') == -1:
                    directory_list.append(iterable[0])
    return directory_list


def get_python_files(directory_list):
    python_file_paths = []
    for directory in directory_list:
        files = os.listdir(directory)
        for file in files:
            if file.endswith('.py'):
                python_file_paths.append(os.path.join(directory, file))
    return python_file_paths


def find_and_replace(python_file_paths):
    for python_file in python_file_paths:
        with open(python_file) as file:
            lines = file.readlines()
            new_lines = []
            for line in lines:
                replacements = []
                if not line.startswith('import') and not line.startswith('from') and not line.startswith('class'):
                    # couldn't be bothered to adjust for native imports (fixed manually)
                    matches = re.findall(r'\w*', line)
                    new_matches = []
                    if len(matches) != 0:
                        for match in matches:
                            if match != '' and len(re.findall(r'[a-z][A-Z]', match)) != 0:
                                new_matches.append(match)
                        if len(new_matches) != 0:
                            for match in new_matches:
                                new_name = camel_to_snake(match)
                                replacements.append([match, new_name])
                new_line = line
                if len(replacements) != 0:
                    for row in replacements:
                        new_line = new_line.replace(row[0], row[1])
                new_lines.append(new_line)
            with open(python_file, 'w') as new_file:
                new_file.writelines(new_lines)


def camel_to_snake(name):
    name = re.sub('(.)([A-Z][a-z]+)', r'\1_\2', name)
    return re.sub('([a-z0-9])([A-Z])', r'\1_\2', name).lower()


# execute code
repository = '/Users/amcg0011/GitRepos/msbrainpy/msbrainpy'
relevant_directories = get_relevant_directories(repository)
python_files = get_python_files(relevant_directories)
find_and_replace(python_files)
