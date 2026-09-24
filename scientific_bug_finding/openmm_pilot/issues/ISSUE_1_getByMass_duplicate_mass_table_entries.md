# `Element.getByMass` can never return elements whose tabulated mass exactly matches another element's

## Environment

- OpenMM 8.6 (`openmm/app/element.py`)
- Python 3.11

## Description

`Element.getByMass` is documented to return "the element whose mass is
CLOSEST to the requested mass," but for several elements this is never
possible, regardless of the input, because the internal cache used by
`getByMass` is a plain `dict` keyed by the raw mass value:

```python
if Element._elements_by_mass is None:
    Element._elements_by_mass = OrderedDict()
    for elem in sorted(Element._elements_by_symbol.values(), key=lambda x: x.mass):
        Element._elements_by_mass[elem.mass.value_in_unit(daltons)] = elem
```

Several pairs of elements in `element.py`'s built-in table share the exact
same tabulated mass (both are unstable synthetic elements whose masses are
given as rounded whole-number estimates rather than measured values):

- berkelium and curium: both `247` daltons
- dubnium and lawrencium: both `262` daltons

When the cache dict is built, the second entry for a given mass silently
overwrites the first, so one of each pair is dropped from
`_elements_by_mass` entirely. `getByMass` can therefore never return
curium or lawrencium, no matter how close the requested mass is to their
tabulated value -- it always returns berkelium or dubnium instead, even
for input values strictly closer to curium/lawrencium than to
berkelium/dubnium's own tabulated mass.

## Reproduction

```python
>>> from openmm.app import element
>>> from openmm import unit
>>> element.Element.getByMass(247.0 * unit.daltons)
<Element berkelium>
>>> element.Element.getByMass(246.9 * unit.daltons)   # strictly closer to Cm/Bk tie, still Bk
<Element berkelium>
>>> element.Element.getByMass(262.0 * unit.daltons)
<Element dubnium>
>>> element.Element.getByMass(261.8 * unit.daltons)
<Element dubnium>
```

Curium (`Element.getBySymbol('Cm')`) and lawrencium
(`Element.getBySymbol('Lr')`) are perfectly valid, correctly registered
`Element` objects in every other respect -- `getBySymbol` and
`getByAtomicNumber` both work fine for them. Only `getByMass` is affected,
because it is the only lookup path that uses this dict-keyed-by-mass
cache.

## Expected behavior

`getByMass` should return whichever of two same-mass elements is actually
closest to the query, and when two elements are truly tied it should have
some defined (or at least stable/reasonable) tie-breaking behavior --
neither should be permanently unreachable.

## Root cause

`_elements_by_mass` is built as `{mass_value: element}` with no handling
for two elements sharing an identical `mass_value`; the dict construction
silently loses the first entry inserted for any repeated key. A cache
keyed on mass instead needs to store a list of elements per mass (or use
a different data structure entirely, e.g. keeping the existing
mass-sorted list without collapsing it into a dict), since element mass
is not a unique key across the whole table.
