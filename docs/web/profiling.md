---
title: Profiling system in `MGIS`
author: Julien Rigal, Thomas Helfer
date: 2026
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

This section describes the profiling system introduced in `MGIS`, which is meant to:

1. provide a hierarchical and precise measurement of execution times across the library.
2. ensure minimal to zero overhead when disabled, preserving the high-performance computing requirements of the code.

By tying the profiling system directly to the `Context` object, `MGIS` avoids the need to pass an additional profiler object through the call stack. A function using this profiling strategy is recognizable by the presence of a `Context` object and the use of dedicated profiling macros.

# The `ProfilingData` structure

The result of the profiling is not a flat list of timers, but a hierarchical tree. Each node in this tree is represented by the `ProfilingData` structure, which stores:

- `name`: A `std::string` representing the name of the profiled section.
- `time_in_seconds`: A `double` accumulating the total time spent in this section.
- `calls`: An `unsigned int` counting the number of times this section was executed.
- `children`: A `std::vector` of `std::shared_ptr<ProfilingData>`, representing the nested profiled sections called within the current one.

# Managing the Profiling State

The profiling system is disabled by default to guarantee zero overhead in production runs. It can be dynamically controlled through the `Context` object using the following methods:

- `enableProfiling(bool)`: Turns the global profiling on or off for the given context.
- `isProfilingEnabled()`: Returns a boolean indicating the current state of the profiler.
- `getProfilingResultTree()`: Retrieves the root node of the profiling tree (`ProfilingData`) containing all the collected metrics.

# Instrumenting the code

The profiling relies on the **RAII** (Resource Acquisition Is Initialization) idiom. Instead of manually starting and stopping timers, the developer creates a scoped object that starts a timer upon construction and stops it, while updating the tree, upon destruction.

To ensure high readability and ease of use, this mechanism is wrapped in macros.

## The `CatchTimeSection` macro

The `CatchTimeSection` macro is the standard way to profile a block of code. It takes two arguments:
1. The `Context` object.
2. A string literal representing the name of the section.

If profiling is enabled in the provided `Context`, the macro will automatically find its place in the hierarchical tree, start the timer, and record the elapsed time at the end of the current scope. If profiling is disabled, the macro does nothing.

## The `CatchLocalTimeSection` macro

In some specific debugging scenarios, a developer might want to force the profiling of a specific section even if the global profiling state is disabled. The `CatchLocalTimeSection` macro takes a third boolean argument:

~~~~{.cxx}
// The third argument 'true' forces the profiling of this specific scope
CatchLocalTimeSection(ctx, "ForcedSection", true); 
~~~~

# Aggregation and Loops

To prevent the profiling tree from growing indefinitely and consuming too much memory, the system automatically aggregates repeated calls to the same section within the same parent scope. 

If a profiled section is called multiple times (for example, inside a `for` or `while` loop), the profiler does not create a new child node for each iteration. Instead, it finds the existing child node with the same name, increments its `calls` counter, and adds the elapsed time to `time_in_seconds`.

# Example of usage

The following code illustrates how to instrument a function and how the tree hierarchy and loop aggregation behave:

~~~~{.cxx}
void performComputation(mgis::Context& ctx)
{
  ctx.enableProfiling(true);

  // Start a root section
  {
    CatchTimeSection(ctx, "IntegrationStep");
    
    // Nested section
    {
      CatchTimeSection(ctx, "Initialization");
      // ... initialization code ...
    }

    // Loop with a nested section
    for (int i = 0; i < 1000; ++i) {
      CatchTimeSection(ctx, "NewtonRaphsonIteration");
      // ... solver code ...
    }
  }

  // Retrieve and analyze the results
  const auto& root = ctx.getProfilingResultTree();
  // 'root' contains 1 child: "IntegrationStep"
  // "IntegrationStep" contains 2 children: "Initialization" and "NewtonRaphsonIteration"
  // "NewtonRaphsonIteration" will show: calls = 1000
}
~~~~

In this example, despite the `NewtonRaphsonIteration` section being timed 1000 times, only one node is created in the tree under `IntegrationStep`, with its `calls` attribute set to `1000` and its `time_in_seconds` representing the total time spent across all iterations.

> **Note**
> 
> Exhaustive examples and unit tests of the profiling system, including edge cases, can be found in the `tests/ProfilingTest.cxx` file.