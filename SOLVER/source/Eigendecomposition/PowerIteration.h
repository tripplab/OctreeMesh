// PowerIteration.h
// Copyright (C) 2013 Miguel Vargas (miguel.vargas@gmail.com)
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Library General Public
// License as published by the Free Software Foundation; either
// version 2 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public
// License along with this library; if not, write to the Free
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

#ifndef _PowerIteration_h_
#define _PowerIteration_h_

#include <Basic/Float.h>
#include <Basic/Log.h>
#include <Basic/Memory.h>
#include <Container/CSRMatrix.h>
#include <Container/Vector.h>
#include <Solver/Cholesky.h>
#include <math.h>


template <typename T>
int PowerIteration(const CSRMatrix<T>& A, Vector<Vector<T> >& v, T& l, T tolerance, int max_steps, int threads) throw(Memory::Exception)
{
	try
	{
		Log(1, "PowerIteration:");
		Log(1, "-Tolerance:   %.5e", tolerance);
		Log(1, "-MaxSteps:    %i", max_steps);

		Vector<T> w;

		for (int s = 1; s <= max_steps; ++s)
		{
			T ww = 0;
			#pragma omp parallel for default(shared) reduction(+:ww) num_threads(threads)
			for (int i = 1; i <= n; ++i)
			{
				register int* __restrict A_index_i = A.index[i];
				register T* __restrict A_entry_i = A.entry[i];

				register T sum = 0.0;
				int k_max = A.Count(i);
				for (register int k = 1; k <= k_max; ++k)
				{
					sum += A_entry_i[k]*w.entry[A_index_i[k]];
				}
				w.entry[i] = sum;
				ww += sum*sum;
			}

			T norm_ww = sqrt(ww);
			for (int i = 1; i <= n; ++i)
			{
				v.entry[i] = w.entry[i]/norm_ww;
			}
			if (norm_ww < tolerance)
			{
				Log(1, "-Total steps: %i", step);
				break;
			}
		}

		if (step >= max_steps)
		{
			Log(1, "-[Error] PowerIteration did not converge in %i steps", max_steps);
			return -1;
		}

		return step;
	}
	catch (Exception&)
	{
		ReThrow();
	}
}

#endif
