import os, shutil, tempfile

'''A file object to `/dev/null`.'''
DEV_NULL = open(os.devnull, "r+")

def read_file(filename):
    '''Reads a file and returns its contents as a string.'''
    with open(filename) as f:
        return f.read()

def write_file(filename, contents):
    '''Writes a file with the given string as contents.'''
    with open(filename, "w") as f:
        f.write(contents)

class TmpFile:
    '''A context object that represents a temporary file.  The file is deleted
    upon exiting the context.'''

    def __init__(self, contents, suffix=""):
        '''Creates a temporary file with the given string as contents and
        returns it as a `TmpFile`.'''
        self.fd, self.filename = tempfile.mkstemp(suffix=suffix)
        with os.fdopen(self.fd, "w") as f:
            f.write(contents)

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        try:
            os.remove(self.filename)
        except:
            pass

class TmpDir:
    '''A context object that represents a temporary directory.  The file is
    deleted upon exiting the context.'''

    def __init__(self, dir=None):
        '''Creates a temporary directory.'''
        self.dirname = tempfile.mkdtemp(dir=dir)

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        try:
            shutil.rmtree(self.dirname)
        except:
            pass
