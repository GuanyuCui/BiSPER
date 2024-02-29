#ifndef ALGORITHMS_H
#define ALGORITHMS_H

#include "graph.h"

class Algorithms
{
	public:
		using size_type = uint64_t;
		using length_type = uint64_t;

		Algorithms(const Graph & G);
		virtual ~Algorithms() = default;

	protected:
		const Graph & G;
};

#endif