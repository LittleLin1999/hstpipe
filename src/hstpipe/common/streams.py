"""Capture legacy numerical-stage progress without changing its implementation."""
from contextlib import contextmanager,redirect_stdout,redirect_stderr
import sys
class Tee:
    def __init__(self,console,file):self.console,self.file=console,file
    def write(self,text):
        self.console.write(text);self.file.write(text);self.file.flush();return len(text)
    def flush(self):self.console.flush();self.file.flush()
@contextmanager
def stage_log(path):
    with open(path,'a',buffering=1) as file:
        with redirect_stdout(Tee(sys.stdout,file)),redirect_stderr(Tee(sys.stderr,file)):
            yield
