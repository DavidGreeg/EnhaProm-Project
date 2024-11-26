import io, textwrap
from contextlib import redirect_stdout

def outputwrap1(output_func, set_width = 50):
	output = str(output_func)
	wrap_output = textwrap.fill(output,
                                width = set_width)
	print(wrap_output)

def outputwrap(output_func, width = 50, 
               args = (), kwargs = None):
    if kwargs is None:
        kwargs = {}
    output_stream = io.StringIO()
    with redirect_stdout(output_stream):
        func_df = output_func(*args, **kwargs)
    full_output = output_stream.getvalue()
    wrapped_output = textwrap.fill(full_output,
                                   width = width)
    subst_output = wrapped_output.replace(". ", ".\n\n")
    print(subst_output)
    print(func_df)
