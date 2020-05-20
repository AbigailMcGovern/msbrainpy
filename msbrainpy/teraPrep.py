import os
import re
import shutil
from collections import defaultdict
from pathlib import Path


# --------------------------------- Correctly place and name files for tera_stitcher  -----------------------------------

def write_tera_files(source_dir, out_dir, prefix, left_ls='ultra_ii Filter0000.tif',
                   right_ls='ultra_ii Filter0001.tif', channels=r"_C[0-9]{2}_xyz", chan_name=r"C[0-9]{2}",
                   left_tiles=(" x 00", " x 01"), right_tiles=(" x 02", " x 03"), x_tiles=4, y_tiles=4,
                   X=2.032, Y=2.032, Z=3, x_width=2160, y_length=2560, overlap=.2,
                   tile_string="[0{0} x 0{1}]", z_slice=r"xyz-Table\s_z\d{4}"):
    source_dir = Path(source_dir)
    out_dir = Path(out_dir)
    new_name = prefix + '_TS'
    new_path = out_dir / new_name
    if not new_path.exists():
        os.mkdir(new_path)
    x_size = X * x_width
    y_size = Y * y_length
    inv_ol = 1 - overlap
    filelist = os.listdir(source_dir)
    lightsheets = get_correct_files(filelist, left_ls, right_ls, left_tiles, right_tiles)
    channels_dict = get_channels(lightsheets, new_path, channels, chan_name)
    write_channel_files(channels_dict, source_dir, x_tiles, y_tiles, y_size, x_size, Z, inv_ol, z_slice, tile_string)


def get_correct_files(filelist, left_ls, right_ls, left_tiles, right_tiles):
    left_lightsheet = []
    right_lightsheet = []
    for _file in filelist:
        if _file.endswith(left_ls):
            left_lightsheet.append(_file)
        if _file.endswith(right_ls):
            right_lightsheet.append(_file)
    files = []
    l = 0
    for tile in left_tiles:
        for _file in left_lightsheet:
            if _file.find(tile) != -1:
                files.append(_file)
                l += 1
    r = 0
    for tile in right_tiles:
        for _file in right_lightsheet:
            if _file.find(tile) != -1:
                files.append(_file)
                r += 1
    print('{} and {} files were found for the left and right lightsheets, respectively'.format(l, r))
    return files


def get_channels(lightsheets, new_path, channels, chan_name):
    unique_channels = []
    for _file in lightsheets[0]:
        a = re.search(channels, _file)
        if a.group(0) not in unique_channels:
            unique_channels.append(a.group(0))
    print('Unique channels found: {}'.format(len(unique_channels)))
    new_name = new_path.stem
    channel_dict = defaultdict()
    channel_paths = []
    for channel in unique_channels:
        b = re.search(chan_name, channel)
        b = b.group(0)
        channel_name = new_name + "_" + b
        channel_path = new_path / channel_name
        channel_paths.append(channel_path)
        os.mkdir(channel_path)
        channel_files = []
        for _file in lightsheets:
            if _file.find(channel) != -1:
                channel_files.append(_file)
        print('{} files were found for channel {}'.format(len(channel_files), b))
        print('Only appropriately sided tiles were added')
        channel_dict[b] = defaultdict()
        channel_dict[b]['chan_string'] = channel
        channel_dict[b]['files'] = channel_files
        channel_dict[b]['path'] = channel_path
    return channel_dict


def write_channel_files(channels_dict, source_dir, x_tiles, y_tiles, y_size, x_size, Z, inv_ol, z_slice, tile_string):
    for channel in channels_dict.keys():
        channel_files = channels_dict[channel]['files']
        channel_path = channels_dict[channel]['path']
        for row in range(y_tiles):
            micron10th_position = y_size * inv_ol * row * 10
            micron10th_position = round(micron10th_position)
            micron10th_position = "{:06d}".format(micron10th_position)
            new_path = channel_path / micron10th_position
            os.mkdir(new_path)
        rows_in_channel = os.listdir(channel_path)
        for row in rows_in_channel:
            row_dir = channel_path / row
            for tile in range(x_tiles):
                micron10th_position = x_size * inv_ol * tile * 10
                micron10th_position = round(micron10th_position)
                micron10th_position = "{:06d}".format(micron10th_position)
                new_dir = row + "_" + micron10th_position
                new_path = row_dir / new_dir
                os.mkdir(new_path)
        for i in range(y_tiles):
            row_dir = channel_path / rows_in_channel[i]
            cols_in_row = os.listdir(row_dir)
            for j in range(x_tiles):
                col_dir = row_dir / cols_in_row[j]
                tile = tile_string.format(i, j)
                files_to_copy = []
                for _file in channel_files:
                    if _file.find(tile) != -1:
                        files_to_copy.append(_file)
                for _file in files_to_copy:
                    file_source = source_dir / _file
                    file_dest = col_dir / _file
                    shutil.copy2(file_source, file_dest)
                stack_files = os.listdir(col_dir)
                for Zslice in stack_files:
                    a = re.search(z_slice, Zslice)
                    slice_no = re.search(r"\d{4}", a.group(0))
                    slice_no = int(slice_no.group(0))
                    position = round(slice_no * Z * 10)
                    pos_string = '{:06d}'.format(position)
                    new_name = pos_string + '.tif'
                    old_adress = col_dir / Zslice
                    new_adress = col_dir / new_name
                    os.rename(old_adress, new_adress)
                print('{} files were saved to {}. \n_files were renamed according to position'.format(len(stack_files),
                                                                                                     col_dir))


# ---------------------- Correctly place and name files + remove prior files for tera_stitcher  -------------------------
# Somewhere in the middle of refactoring the write_terra_files() function, never finished this.
#    Will return when (and if) inclined.
def make_tera_dirs(channels_dict, y_tiles, y_size, inv_ol):
    for channel in channels_dict.keys():
        channel_files = channels_dict[channel]['files']
        channel_path = channels_dict[channel]['path']
        for row in range(y_tiles):
            micron10th_position = y_size * inv_ol * row * 10
            micron10th_position = round(micron10th_position)
            micron10th_position = "{:06d}".format(micron10th_position)
            new_path = channel_path / micron10th_position
            os.mkdir(new_path)
        rows_in_channel = os.listdir(channel_path)
        for row in rows_in_channel:
            row_dir = channel_path / row
            for tile in range(x_tiles):
                micron10th_position = x_size * inv_ol * tile * 10
                micron10th_position = round(micron10th_position)
                micron10th_position = "{:06d}".format(micron10th_position)
                new_dir = row + micron10th_position
                new_path = row_dir / new_dir
                os.mkdir(new_path)
    return rows_in_channel


def copy2_tera_dir(rows_in_channel, y_tiles, x_tiles, tile_string):
    for i in range(y_tiles):
        row_dir = channel_path / rows_in_channel[i]
        cols_in_row = os.listdir(row_dir)
        for j in range(x_tiles):
            col_dir = row_dir / cols_in_row[j]
            tile = tile_string.format(i, j)
            files_to_copy = []
            for _file in channel_files:
                if _file.find(tile) != -1:
                    files_to_copy.append(_file)
            for _file in files_to_copy:
                file_source = source_dir / _file
                file_dest = col_dir / _file
                shutil.copy2(file_source, file_dest)
                yield col_dir


def rename_in_stack_dir(col_dir, z_slice, Z):
    stack_files = os.listdir(col_dir)
    for Zslice in stack_files:
        a = re.search(z_slice, Zslice)
        slice_no = re.search(r"\d{4}", a.group(0))
        slice_no = int(slice_no.group(0))
        position = round(slice_no * Z * 10)
        pos_string = '{:06d}'.format(position)
        new_name = pos_string + '.tif'
        old_adress = col_dir / Zslice
        new_adress = col_dir / new_name
        os.rename(old_adress, new_adress)
