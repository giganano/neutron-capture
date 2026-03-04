#!/usr/bin/env python


def binary_ns(t):
	if t:
		return t**-1.1
	else:
		return 0


def tertiary_ns(t, minimum = 3):
	if t >= minimum:
		return t**-1.1
	else:
		return 0


