---
layout: post
title: Using mock to write unit tests for random functions
---

Writing unit tests is important,
but this becomes difficult for functions which use random elements.
One possible approach would be to use known seeds for the random number generation,
allowing you to have a deterministic sequence of random numbers,
so that regression tests could be written.
Recently I've been writing unit tests for a function which uses a sequence of `random.random()` calls to determine what is returned
(this is for mutating individuals within a genetic algorithm).
Because of this, I want full control over what is returned by `random.random()` so that I can test all possible permutations of what happens inside the function.

We can use the [mock package](https://docs.python.org/3/library/unittest.mock.html) to replace the behaviour of any calls to `random.random()` with behaviour we define ourselves.
For my purposes, I wanted to make the calls to random be replaced with a sequence of predefined values.
This can be done through passing a generator to the `side_effect` keyword argument of a `mock.Mock` object.
This small test shows how calls to `random.random` are replaced by the list of numbers I provided.

```python
import mock  # in python 3: `from unittest import mock` will work
import random


def RandMock(randseq):
	return mock.Mock(side_effect=(val for val in randseq))

def test_randomness():
	with mock.patch('random.random', RandMock([1, 2, 3, 4])):
		assert random.random() == 1
		assert random.random() == 2
		assert random.random() == 3
		assert random.random() == 4
```


Moving on to a more real example,
I had a function which performed a [blend crossover](http://www.tomaszgwiazda.com/blendX.htm) between parents
to produce offspring with values which are a mix of the two parents (or slightly outside).
This had a call to `random.random()` to first decide if a crossover is performed (90% chance),
and then if so, two more calls (one per child) to decide where in the range of values between the two parents the child value is
(1.0 would extend slightly beyond the highest value, 0.5 would choose a value directly between the parents' values).
Here I've also used a [`pytest.mark.parametrize`](https://docs.pytest.org/en/latest/parametrize.html#pytest-mark-parametrize-parametrizing-test-functions)
decorator to iterate over many different test cases with a single function.


```python
import mock
import random
import pytest


## Somewhere in my package....
bounds = ((10.0, 30.0),)
blend_alpha = 0.1
blend_probability = 0.9

def clamp(val, lower, upper):
	"""He's champin' for a clampin'!"""
	val = max(val, lower)
	val = min(val, upper)
	return val


def blend_crossover(candidates):
	new_candidates = []

	for mother, father in zip(candidates[::2], candidates[1::2]):
		brother, sister = [], []
		for a, b, (min_bound, max_bound) in zip(mother, father, bounds):
			if random.random() < blend_probability:
				smallest, largest = min(a, b), max(a, b)
				# range between parents
				width = largest - smallest
				# amount we go out of bounds from natural range
				extra = width * blend_alpha

				# start point is (smallest - delta)
				# we then move up to (width + 2 * extra) from this point
				a = ((smallest - extra) +
					 random.random() * (width + 2 * extra))
				b = ((smallest - extra) +
					 random.random() * (width + 2 * extra))
				a = clamp(a, min_bound, max_bound)
				b = clamp(b, min_bound, max_bound)

			brother.append(a)
			sister.append(b)

		new_candidates.append(tuple(brother))
		new_candidates.append(tuple(sister))

	return new_candidates


## In my test suite....
def RandMock(randseq):
	return mock.Mock(side_effect=(val for val in randseq))


@pytest.mark.parametrize('randseq,mum,dad,child1,child2', [
	# randseq - (roll for performing crossover,
	#            child1 position, child2 position)
	((0.0, 0.0, 0.0), (15.0,), (25.0,), (14.0,), (14.0,)),  # check minimum crossover
	((0.0, 0.0, 1.0), (15.0,), (25.0,), (14.0,), (26.0,)),  # check other limit
	((0.0, 0.5, 0.5), (15.0,), (25.0,), (20.0,), (20.0,)),  # check no crossover
	((1.0, 1.0, 1.0), (15.0,), (25.0,), (15.0,), (25.0,)),  # check no crossover
])
def test_crossover(randseq, mum, dad, child1, child2):
	with mock.patch('random.random', RandMock(randseq)):
		ret = blend_crossover([mum, dad])

		assert ret[0] == child1
		assert ret[1] == child2

```
