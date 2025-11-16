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
`MGIS` Version 3.1.

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

1.  a -string
2.  a -string
3.  a couple of a function allowing to return an error message from an
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
`std::source_location::current_location`
method, is automatically added to the error message returned by the
`getErrorMessage` method.

Note that this feature is currently only supported by `gcc` compilers.

In some cases, the source location is not meaningful. In this case, the
`registerErrorMessageWithoutSourceLocation` method can be used.

If the information about the source location are not wanted, the
`getRawErrorMessage` method can be used.

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

-   exceptions
-   having a boolean data member stating if the object is valid (this
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

# Interaction with an external library that may use exceptions

## The `registerExceptionInErrorBacktrace` function

Call to external libraries that relies on the usage of exceptions must
be encapsulated in appropriate `try/catch`
blocks as follows:

``` cpp
try{
  ....
} catch(...)
  registerExceptionInErrorBacktrace(e);
}
```

The `registerExceptionInErrorBacktrace` is
a Lippincott-like helper function which translate exceptions derived
from `std::exception` into error messages.

If the external library to be used, used another exception hierarchy,
then appropriate versions of this helper function shall be created.
