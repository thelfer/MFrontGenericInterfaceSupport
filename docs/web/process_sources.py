# -*- coding: utf-8 -*-


import sys
import os
import shutil
import zipfile


def process(format="md"):


    # Directories to scan
    demos = "../../demos/mgis_fenics"
    retval = os.getcwd()

    # Create directory in documentation tree for demo
    doc_demo_dir = './'
    if not os.path.exists(doc_demo_dir):
        os.makedirs(doc_demo_dir)

    supp_file_path = os.path.join(demos, "./supplementary_files/")
    copied_supp_file_path = os.path.join(doc_demo_dir, "./supplementary_files/")
    if not os.path.exists(copied_supp_file_path):
        os.makedirs(copied_supp_file_path)
    tocopy_files = [f for f in os.listdir(supp_file_path)]
    print("Files to copy:", tocopy_files)

    for f in tocopy_files:
        source = os.path.join(supp_file_path, f)
        shutil.copy(source, copied_supp_file_path)

    ipynb_files =  [f for f in os.listdir(demos) if f.endswith("ipynb") ]
    for f in ipynb_files:
        ff = os.path.join(demos, f)
        shutil.copy(ff, doc_demo_dir)
        output = os.path.join(doc_demo_dir, os.path.splitext(f)[0])
        print(ff)
        print(output)

    os.chdir(doc_demo_dir)
    for f in ipynb_files:
        # # try first with latex environments
        # ret = os.system("jupyter-nbconvert --to html_with_lenvs %s --output %s" % (ff,  os.path.join("../doc/", output)))
        # if not ret == 0:
        try:
            # otherwise export without latex environments
            if format == "html":
                ret = os.system("jupyter-nbconvert --to html %s --output %s" % (ff,  os.path.join("../doc/", output)))
            elif format == "md":
                fmd = f.replace(".ipynb", ".md")
                fhtml = f.replace(".ipynb",".html")
                ret = os.system("jupyter-nbconvert --to markdown {} --output {}".format(f, fmd))
                change_ipynb_links_to_html(fmd)
                replace_headers(fmd)
                # print(fmd, fhtml)
                #ret = os.system("pandoc -s --css ../main.css -f markdown-markdown_in_html_blocks --mathjax {} -o {}".format(fmd, fhtml))
        except:
            raise RuntimeError("Unable to convert ipynb file to a .py ({})".format(f))

        # remove_ipython_magic(output+".py")

def remove_ipython_magic(file):
    with open(file, "r+") as f:
        new_f = f.readlines()
        f.seek(0)
        for line in new_f:
            if "get_ipython().run_line_magic" not in line:
                f.write(line)
        f.truncate()


from tempfile import mkstemp
from shutil import move
from os import remove

def replace(source_file_path, pattern, substring, startwith=False):
    fh, target_file_path = mkstemp()
    with open(target_file_path, 'w') as target_file:
        with open(source_file_path, 'r') as source_file:
            for line in source_file:
                if type(pattern) == str:
                    target_file.write(line.replace(pattern, substring))
                else:
                    for (p, s) in zip(pattern, substring):
                        if startwith:
                            b = line.startswith(p)
                        else:
                            b = p in line
                        if b:
                            target_file.write(line.replace(p, s))
                            break
                    else:
                        target_file.write(line)

    remove(source_file_path)
    move(target_file_path, source_file_path)

def change_ipynb_links_to_html(file):
    replace(file, ".ipynb)", ".html)")

def replace_headers(file):
    replace(file, ("# ", "## ", "### "), ("% ", "# ", "## "), startwith=True)

if __name__ == '__main__':
    process()