---
title: em-widgets demo
---

These are the same notebooks as in `test/index.html`, but mounted through MyST's
`{anywidget}` directive — i.e. exactly how they appear on a real site.

## SEM ray diagram

:::{anywidget} ../observable-notebook.mjs
{
  // The body is JSON5: comments and unquoted keys are fine.
  notebook: "/notebooks/sem-ray-diagram.html",
  cells: ["rayDiagram"]   // omit to render the whole notebook, prose and all
}
:::

## Geometric aberrations

:::{anywidget} ../observable-notebook.mjs
{
  notebook: "/notebooks/geometric-aberrations.html"
}
:::

## Aperture autocorrelation

:::{anywidget} ../observable-notebook.mjs
{
  notebook: "/notebooks/aperture-autocorrelation.html"
}
:::

## Probe aberrations

:::{anywidget} ../observable-notebook.mjs
{
  notebook: "/notebooks/probe-aberrations.html"
}
:::

## STEM measurements (data loaded from a zarr ZipStore)

:::{anywidget} ../observable-notebook.mjs
{
  notebook: "/notebooks/stem-measurements.html"
}
:::
