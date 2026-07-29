---
title: Error handling in `MGIS`
author: Thomas Helfer
date: 2025
lang: en-EN
numbersections: true
documentclass: article
from: markdown+tex_math_single_backslash
geometry:
  - margin=2cm
papersize: a4
link-citations: true
colorlinks: true
figPrefixTemplate: "$$i$$"
tabPrefixTemplate: "$$i$$"
secPrefixTemplate: "$$i$$"
eqnPrefixTemplate: "($$i$$)"
bibliography: bibliography.bib
---

This section describes the error handling strategy introduced in
`MGIS` Version 3.1 which is meant to:

1.  ensure performances and compatibility with the HPC techniques used
    in the code (MPI currently and eventually OpenMP in the future).
2.  provide rich error messages to the end-user.

The first point eliminates the use of exceptions in most cases (see
below the special cases of constructors).

Some other standard strategies provided by the standard have also been
discarded, although their design have inspired the proposed solution.
For example, `std::error_code` only allows
to use predefined error messages that are not able to describe the
context.

A function or a method using this error handling strategy is
recognizable as follows:

- It must take a class derived from `AbstractErrorHandler&` as its first
  parameter. Two main classes, derived from `AbstractErrorHandler&`, are
  available: `ContractViolationHandler` and `Context`. Consequently, if
  a function has a `AbstractErrorHandler&` as its first parameter, the
  reader knows that it may fail, and that the return type of the
  function should be interpreted according to the rule below.
- It returns a value that can be invalid. The most common return type is
  a boolean which indicates success (`true`) or failure (`false`).
- It must be declared `noexcept`.

# Returned value of functions (or methods) that may fail

A function that may fail shall generally return either:

- a boolean associated with the success of the function, i.e. if a
  function returns `true`, it succed.
- an `std::optional` object holding the results if those results are
  stored on the stack.
- an `std::unique_ptr` or an `std::shared_ptr` object holding if those
  results are stored on the heap. In this case, the failure of the
  function is indicated by the fact that the underlying pointer is null.
- a class to which an `InvalidValue` object is implicitly convertible.

> ** Note **
> 
> Those returned value shall be mandatory marked with the `[[nodiscard]]`
> attribute. However, this conflicts in many case with the `MGIS_EXPORT`
> attribute. This is a defect of `C++-20`.

# `ContractViolationHandler` and `Context`

Both classses `ContractViolationHandler` and `Context` inherit from
`AbstractErrorHandler`.

## `ContractViolationHandler`

`ContractViolationHandler` is designed to report a contract violation,
i.e. errors that shall not occur by design.

`ContractViolationHandler` can be used in a `constexpr` context, if no
contrat violation is detected. If a contract violation is detected, a
compile-time error is generated since `registerErrorMessage` is not
`constexpr`

### Reporting error to end-user

The `ContractViolationHandler` class provides the `registerErrorMessage`
method which accepts a `C`-string.

The behaviour of this class is to call `std::abort` in case of contract
violation. This behaviour can be changed by passing
`-Denable-exceptions=ON` to `cmake` when compiling `MGIS`.

## `Context`

The `Context` class is used for standard error management.

### Reporting error to end-user

The `Context` class inherits from the `ErrorBacktrace` class, which has
been designed to store error messages in a hierarchical way from the
lowest level of the code up to the highest level function. The
`ErrorBacktrace` is meant to be used through the `Context` object which
is passed as the first argument to most functions.

The `ErrorBacktrace` class, and the `Context` class mostly provides the
`registerErrorMessage` method that can register an error in the form of:

1. a `C`-string
2. a  `C++`-string
3. a couple of a function allowing to return an error message from an
  integer value. This kind of function is provided by many HPC
  libraries.

For convenience, the `registerErrorMessage` always returns an invalid
value, i.e. a value that is convertible to any of the returned type
described in the previous section.

By default, the error messages are packed up to the moment when error(s)
must be reported to the end-user. The error messages can then be
retrieved by the `getErrorMessage` method (or `getRawErrorMessage`, see
below). However, if `MGIS` is compiled with the flag
`-Denable-exceptions=ON`, an exception is thrown instead using
`mgis::raise`.

### Source location

In debug mode, the source location, as returned by the
`std::source_location::current_location` method, is automatically added
to the error message returned by the `getErrorMessage` method.

Note that this feature is currently only supported by `gcc` compilers.

In some cases, the source location is not meaningful. In this case, the
`registerErrorMessageWithoutSourceLocation` method can be used.

If the information about the source location are not wanted, the
`getRawErrorMessage` method can be used.

## Example of usage of the `registerErrorMessage` method

The following code illustrates the usage of the `ErrorBacktrace`,
through a `Context` object:

~~~~{.cxx}
void processNext(){};

bool f3(Context &ctx)
{
  // for this example, the message is useless
  // In pratice, one shall report the cause of the error
  // and not expose details, like the function name.
  //
  // Examples:
  //
  // - "negative temperature detected"
  // - "non convergence of the nonlinear solver"
  //
  // Note that the function name, the source file and the
  // line number are automatically added in debug mode.
  return ctx.registerErrorMessage("invalid call");
}

bool f2(Context &ctx)
{
  if (!f3(ctx)) {
    // f2 fails, but we don't have any more information
    // to add for the end user (i.e. f3 is an internal 
    // method and a message like `f3 failed` is not 
    // meaningful), so we just return
    return false;
  }
  processNext();
  return true;
}

bool f1(Context &ctx)
{
  if (!f2(ctx)) {
    return e.registerErrorMessage("invalid call to f2");
  }
  processNext();
  return true;
}
~~~~

If the `f1` function is called, the
following error message is generated in release mode:

~~~~ bash
invalid call to f2
* invalid call
~~~~

In debug mode, the following message is generated:

~~~~{.bash}
/home/UserDir/tests/core/error_backtrace_handler/error_backtrace_handler_test.cpp:31: in function 'bool f1(mgis::ErrorBacktrace&)': invalid call to f2
* /home/UserDir/tests/core/error_backtrace_handler/error_backtrace_handler_test.cpp:14: in function 'bool f3(mgis::ErrorBacktrace&)': invalid call
~~~~

### Warning about usage of `C`-strings

Note that in the case of a `C`-string, the string is not copied and only
the pointer is stored. The developer must then ensure that this string
is not destroyed. As a rule of thumb, this string shall belong to the
data section of the binary.

### Warning about usage of `C++`-strings

`C++`-strings are the best way to report context sensitive error message.
However, to reduce code bloat, building a complex error message shall
never be implemented in a template function: one shall create a
dedicated non-template function implemented in a source file.

# The special case of constructors

Constructors don't return values. There are mostly two ways to handle
failure in constructors:

- exceptions
- having a boolean data member stating if the object is valid (this
  strategy is used by the standard `iostream` library for instance).

Here, we propose to use exceptions and to wrap constructors in dedicated
functions.

## The `raise` function

The `raise` function is an utility function to throw exception in a safe
way: building the exception is not done in the `throw` statement. This
is required to avoid a potential undefined behaviour if the constructor
of the exception throws.

The type of the exception thrown is given by the first template argument
of the `raise` function and defaults to `std::runtime_exception`.

## The `construct` function

The `construct` function calls the constructor of an object that may
throw an exception and returns a `std::optional` object that holds the
object if the call to the constructor did not throw an exception.

If the constructor threw an exception, then the error message hold by
the exception is registered in the instance of the `ErrorBacktrace`
class which is passed as the first argument of the function and an empty
optional object is returned.

Aside from the first argument (a reference to an instance of the
`ErrorBacktrace` class), all the other arguments are forwarded to the
constructor of the object.

### The `MGIS_CONSTRUCT` macro

The `MGIS_CONSTRUCT` macro is a wrapper
around the `construct` function which adds
the current source location when required (typically in the `debug`
mode).

### The `MGIS_TRY_CONSTRUCT` macro

A typical pattern of usage of the `construct` function (through the
`MGIS_CONSTRUCT` macro) is to try to build an object and:

1. to defer the result in case of success.
2. to stop the execution of the current function in case of failure.

Here is a typical example of this pattern:

~~~~{.cxx}
auto tmp_v = MGIS_CONSTRUCT(ObjectType, e, ...);
if(!tmp_v.has_value()){
  return false;
}
auto& v = *(tmp_v);
~~~~

The `MGIS_TRY_CONSTRUCT` macro reduces
this code as follows:

~~~~{.cxx}
MGIS_TRY_CONSTRUCT(ObjectType, v, e, ...);
// Here you can work with variable v which is a reference
// to the ObjectType built by the "construct" function
~~~~

## The `make_unique` function

The `make_unique` function tries to allocate an object on the heap and
stores it in a `std::unique_ptr`.

If an exception is thrown during the construction of the object, the
error message held by the exception is registered in an instance of the
`ErrorBacktrace` and an empty pointer is
returned.

### Helper macros for the `make_unique` function

The `MGIS_MAKE_UNIQUE` and
`MGIS_TRY_MAKE_UNIQUE` macros are similar
to the `MGIS_CONSTRUCT` and
`MGIS_TRY_CONSTRUCT` macros respectively.

## The `make_unique_as` function

The `make_unique_as` function is similar to the `make_unique` function
except that the built object is stored in a `std::unique_ptr` of some
base class. This method is useful in a polymorphic context.

### Helper macros for the `make_unique_as` function

The `MGIS_MAKE_UNIQUE_AS` and `MGIS_TRY_MAKE_UNIQUE_AS` macros are
similar to the `MGIS_CONSTRUCT` and `MGIS_TRY_CONSTRUCT` macros
respectively.

## The `make_shared` function

The `make_shared` function is similar to the `make_unique` function
except that the result is stored in a `std::shared_ptr`.

### Helper macros for the `make_shared` function

The `MGIS_MAKE_SHARED` and `MGIS_TRY_MAKE_SHARED` macros are similar to
the `MGIS_CONSTRUCT` and `MGIS_TRY_CONSTRUCT` macros respectively.

## The `make_shared_as` function

The `make_shared_as` function is similar to the `make_unique_as` function
except that the result is stored in a `std::shared_ptr`.

### Helper macros for the `make_shared_as` function

The `MGIS_MAKE_SHARED_AS` and `MGIS_TRY_MAKE_SHARED_AS` macros are
similar to the `MGIS_CONSTRUCT` and `MGIS_TRY_CONSTRUCT` macros
respectively.

# Interaction with an external library that may use exceptions

## The `registerExceptionInErrorBacktrace` function

Call to external libraries that relies on the usage of exceptions must
be encapsulated in appropriate `try/catch`
blocks as follows:

~~~~{.cxx}
try{
  ....
} catch(...)
  registerExceptionInErrorBacktrace(e);
}
~~~~

The `registerExceptionInErrorBacktrace` is a Lippincott-like helper
function which translate exceptions derived from `std::exception` into
error messages.

If the external library to be used, used another exception hierarchy,
then appropriate versions of this helper function shall be created.
