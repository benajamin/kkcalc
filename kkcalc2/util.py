"""Module to provide utility functions and definitions for the kkcalc package."""

from typing import Callable, Any, TypeAlias, ParamSpec, TypeVar

# Setup typing for the decorator function.
P = ParamSpec("P")
T = TypeVar("T")
WrappedFunctionDecorator: TypeAlias = Callable[[Callable[P, T]], Callable[P, T]]


# Decorator for copying docstrings
def doc_copy(
    template_func: Callable[..., Any] | property,
) -> WrappedFunctionDecorator[P, T]:
    """
    Copy the doc string of the given function to the decorated function.

    If a property is passed, the signature of the getter is used by default.
    To use the setter / deleter, use `@doc_copy(prop.setter)`.

    Parameters
    ----------
    template_func : Callable
        Function whose docstring is to be copied.

    Returns
    -------
    WrappedFunctionDecorator
        A decorator that copies the docstring from `template_func` to the decorated function.
    """

    def decorator(f: Callable[P, T]) -> Callable[P, T]:
        """
        Apply the docstring from the decorator parameter `template_func`.

        Parameters
        ----------
        f : Callable
            The function to which the docstring will be copied.

        Returns
        -------
        Callable
            The function `f` with the docstring copied from `template_func`.
        """
        # f is the decorated function
        if (
            isinstance(template_func, property)
            and hasattr(template_func, "fget")
            and hasattr(template_func.fget, "__code__")
        ):
            # Instead use getter method
            tfunc = template_func.fget
        elif hasattr(template_func, "__code__"):
            tfunc = template_func
        else:
            raise ValueError(
                "Template function is not a property or function with a '__code__' attribute."
            )

        template_pnames = tfunc.__code__.co_varnames
        f_pnames = f.__code__.co_varnames

        # Check the function signature matches
        for i, f_pname in enumerate(f_pnames):
            if template_pnames[i] != f_pname:
                raise ValueError(
                    f"Function signature of {f} does not match that of {tfunc}'s parameter '{template_pnames[i]}'."
                )
        # Copy the documentation
        f.__doc__ = tfunc.__doc__
        return f

    return decorator
