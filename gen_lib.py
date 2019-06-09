import shutil,os

def make_dir(path,clear=False,php=False):
    prep_output({0:path},clear=clear,php=php)

def clear_dir(path,clear=True,php=False):
    prep_output({0:path},clear=clear,php=php)

def prep_output(output_dirs={0:'output'},clear=False,width_100=False,img_extra='',php=False):
    if width_100:
        img_extra = "width='100%'"

    txt = []
    txt.append('<?php')
    txt.append('foreach (glob("*.png") as $filename) {')
    txt.append('    echo "<img src=\'$filename\' {img_extra}> ";'.format(img_extra=img_extra))
    txt.append('}')
    txt.append('?>')
    show_all_txt = '\n'.join(txt)

    txt = []
    txt.append('<?php')
    txt.append('foreach (glob("*.png") as $filename) {')
    txt.append('    echo "<img src=\'$filename\' {img_extra}> <br />";'.format(img_extra=img_extra))
    txt.append('}')
    txt.append('?>')
    show_all_txt_breaks = '\n'.join(txt)

    for value in output_dirs.values():
        if clear:
            try:
#                shutil.rmtree(value)
                os.system('rm -rf {}/*'.format(value))
            except:
                pass
        try:
            os.makedirs(value)
        except:
            pass
        if php:
            with open(os.path.join(value,'0000-show_all.php'),'w') as file_obj:
                file_obj.write(show_all_txt)
            with open(os.path.join(value,'0000-show_all_breaks.php'),'w') as file_obj:
                file_obj.write(show_all_txt_breaks)

def adjust_axes(ax_0,ax_1):
    """
    Force geospace environment axes to line up with histogram
    axes even though it doesn't have a color bar.
    """
    ax_0_pos    = list(ax_0.get_position().bounds)
    ax_1_pos    = list(ax_1.get_position().bounds)
    ax_0_pos[2] = ax_1_pos[2]
    ax_0.set_position(ax_0_pos)
