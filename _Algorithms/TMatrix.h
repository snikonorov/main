#ifndef TMATRIX_H_INCLUDED
#define TMATRIX_H_INCLUDED

#include <algorithm>
#include <vector>
#include <cmath>

///----------------------------------------------------------------

struct IteratorNotIterable
{
};

struct IteratorNotDereferencable
{
};

struct IteratorNotIncrementable
{
};

struct IteratorNotDecrementable
{
};

struct IncompatibleIterators
{
};

struct IncompatibleDimensions
{
};

struct InvalidScalarCast
{
};

struct NonSquareMatrix
{
};

struct SingularMatrix
{
};

///----------------------------------------------------------------

// basic memory iterator
template<class Q, class TContainer, class CIterator>
struct TIterator
{
	TContainer& _container;

	const Q* _beg;
	const Q* _end;
	Q* _cur;
	unsigned _step;		//base pointer increment

	TIterator(TContainer& _container, const Q* _beg, const Q* _end, Q* _cur, unsigned _step = 1)
			  : _container(_container), _beg(_beg), _end(_end), _cur(_cur), _step(_step)
	{
	}

	bool operator==(const TIterator& i) const
	{
		return _cur == i._cur;
	}

	bool operator!=(const TIterator& i) const
	{
		return _cur != i._cur;
	}

	void operator++()
	{
		if (_cur < _end)
		{
			_cur += _step;

			if (_cur > _end)
				_cur = (Q*) _end;
		}
		else
			throw IteratorNotIterable();
	}

	void operator++(int)
	{
		if (_cur < _end)
		{
			_cur += _step;

			if (_cur > _end)
				_cur = (Q*) _end;
		}
		else
			throw IteratorNotIterable();
	}

	void operator--()
	{
		if (_cur >= _beg)
		{
			_cur -= _step;

			if (_cur < _beg)
				_cur = (Q*) _beg;
		}
		else
			throw IteratorNotIterable();
	}

	void operator--(int)
	{
		if (_cur >= _beg)
		{
			_cur -= _step;

			if (_cur < _beg)
				_cur = (Q*) _beg;
		}
		else
			throw IteratorNotIterable();
	}

	CIterator operator+(int p) const
	{
		if (_cur != _end)
			return CIterator(
								_container,
								_beg,
								_end,
								(Q*) (_cur + p*_step < _end ? _cur + p*_step : _end),
								_step
							);
		else
			throw IteratorNotIncrementable();
	}

	CIterator operator-(int p) const
	{
		if (_cur != _beg)
			return CIterator(
								_container,
								_beg,
								_end,
								(Q*) (_cur - p*_step >= _beg ? _cur - p*_step : _beg),
								_step
							);
		else
			throw IteratorNotDecrementable();
	}
};

struct BMatrix
{
	BMatrix()
	{
	}

	virtual ~BMatrix()
	{
	}
};

/// TMatrix ~ [internally] array of rows of values
template<class T>
class TMatrix : public BMatrix
{
	//friend struct pmatrix;

	//protected:
	public:

	enum StrFormat
	{
		Plain = 0,
		List = 1
	};

	using type = T;

    /// "pseudo-vector" object
    template<class Q>
    struct pvector
    {
        template<class TContainer>
		struct iterator: TIterator<Q, TContainer, iterator<TContainer>>
		{
			template<class... Args>
			iterator(Args&&... args) : TIterator<Q, TContainer, iterator<TContainer>>(std::forward<Args>(args)...)
			{
			}

			auto& operator*()
			{
				return *this->_cur;
			}

			auto& operator*() const
			{
				return *this->_cur;
			}
		};

		const Q* _beg;
		const Q* _end;
        unsigned _step;

        pvector(const Q* _beg, const Q* _end, unsigned _step)
                : _beg(_beg), _end(_end), _step(_step)
        {
        }

        decltype(auto) begin() const
        {
            return iterator<const pvector<Q>>(*this, _beg, _end, (Q*) _beg, _step);
        }

        decltype(auto) end() const
        {
            return iterator<const pvector<Q>>(*this, _beg, _end, (Q*) _end, _step);
        }

		decltype(auto) begin()
		{
			return iterator<pvector<Q>>(*this, _beg, _end, (Q*)_beg, _step);
		}

		decltype(auto) end()
		{
			return iterator<pvector<Q>>(*this, _beg, _end, (Q*)_end, _step);
		}

		auto& operator[](int idx)
		{
			return *(begin() + idx);
		}

		const auto& operator[](int idx) const
		{
			return *(begin() + idx);
		}

		Q operator&(const pvector<Q>& v) const
		{
			Q res = 0;

			auto ql = begin();
			auto qr = v.begin();
			while(ql != end() && qr != v.end())
			{
				res += (*ql)*(*qr);

				ql++;
				qr++;
			}

			if (ql != end() || qr != v.end())
				throw IncompatibleDimensions();

			return res;
		}

		void operator+=(const pvector<Q>& v) const
		{
			auto ql = begin();
			auto qr = v.begin();
			while (ql != end() && qr != v.end())
			{
				*ql += *qr;

				ql++;
				qr++;
			}

			if (ql != end() || qr != v.end())
				throw IncompatibleDimensions();
		}

		operator std::vector<Q>() const
		{
			std::vector<Q> res(end() - begin());

			unsigned idx = 0;
			for (const auto& e: *this)
			{
				res[idx++] = e;
			}

			return res;
		}

		void operator=(const std::initializer_list<Q>& l)
		{
			auto ql = begin();
			auto qr = l.begin();

			while (ql != end() && qr != l.end())
			{
				*ql = *qr;

				ql++;
				qr++;
			}

			/*if (ql != end() || qr != l.end())
				throw IncompatibleDimensions();*/
		}

        #if defined(_GLIBCXX_IOSTREAM) || defined(_IOSTREAM_)

        template<class T1>
		friend std::ostream& operator<<(std::ostream& out, const pvector<T1>& m)
        {
            for (const auto& e: m)
            {
                out << e << ' ';
            }

            return out;
        }

        #endif
    };

	/// "pseudo-matrix" object
	template<class Q, class TContainer>  // = TMatrix<Q>>
	class pmatrix
	{
		public:

		using type = Q;

		template<class Container>	// pmatrix<Q, TContainer>
		struct iterator: TIterator<Q, Container, iterator<Container>>
		{
			template<class... Args>
			iterator(Args&&... args) : TIterator<Q, Container, iterator<Container>>(std::forward<Args>(args)...)
			{
			}

			pvector<Q> operator*()
			{
				if (this->_cur != this->_end)
					return pvector<Q>(this->_cur,
									  this->_container._transposed ? this->_container._container._data.data() + this->_container._container._data.size() : this->_cur + this->_step,
									  this->_container._element_inc);
				else
					throw IteratorNotDereferencable();
			}
		};

		TContainer& _container;
		bool _transposed;

		const Q* _beg;
		const Q* _end;
		unsigned _base_inc, _element_inc;		//base pointer and element increment

		pmatrix(TContainer& _container, const Q* _beg, const Q* _end, unsigned _base_inc, unsigned _element_inc, bool _transposed)
				: _container(_container), _transposed(_transposed), _beg(_beg), _end(_end), _base_inc(_base_inc), _element_inc(_element_inc)
		{
		}

		decltype(auto) begin() const
		{
			return iterator<const pmatrix<Q, TContainer>>(*this, _beg, _end, (Q*) _beg, _base_inc);
		}

		decltype(auto) end() const
		{
			return iterator<const pmatrix<Q, TContainer>>(*this, _beg, _end, (Q*) _end, _base_inc);
		}

		decltype(auto) begin()
		{
			return iterator<pmatrix<Q, TContainer>>(*this, _beg, _end, (Q*)_beg, _base_inc);
		}

		decltype(auto) end()
		{
			return iterator<pmatrix<Q, TContainer>>(*this, _beg, _end, (Q*)_end, _base_inc);
		}

		pvector<Q> operator[](int idx)
		{
			return *(begin() + idx);
		}

		decltype(auto) rows()
		{
			return *this;
		}

		decltype(auto) rows() const
		{
			return *this;
		}

		decltype(auto) cols() const
		{
			return _transposed ? _container.rows() : _container.cols();
		}

		std::pair<unsigned, unsigned> Dim() const
		{
			return _transposed ? std::pair<unsigned, unsigned>{_container.NCols(), _container.NRows()} : _container.Dim();
		}

		unsigned NCols() const
		{
			return _transposed ? _container.NRows() : _container.NCols();
		}

		unsigned NRows() const
		{
			return _transposed ? _container.NCols() : _container.NRows();
		}

		protected:

		Q& operator()(unsigned idx)
		{
			return _transposed ? _container._data[idx*_container.NCols() < _container._data.size() ? idx * _container.NCols()
																								   : idx * _container.NCols() % (_container._data.size() - 1)]
							   : _container._data[idx];
		}

		/*const Q& operator()(unsigned idx) const
		{
			return _transposed ? _container._data[idx*_container.NCols() < _container._data.size() ? idx * _container.NCols()
																								   : idx * _container.NCols() % (_container._data.size() - 1)]
							   : _container._data[idx];
		}*/

		public:

		const Q& at(unsigned idx) const
		{
			return _transposed ? _container._data[idx*_container.NCols() < _container._data.size() ? idx * _container.NCols()
																								   : idx * _container.NCols() % (_container._data.size() - 1)]
							   : _container._data[idx];
		}

		operator TMatrix() const
		{
			TMatrix res(Dim());

			/*unsigned idx = 0;
			for (const auto& row : *this)
			{
				for (const auto& e : row)
				{
					res._data[idx++] = e;
				}
			}*/

			for (unsigned idx = 0; idx < res._data.size(); idx++)
			{
				res._data[idx] = this->at(idx);
			}

			return res;
		}

		///---------------------------------------------------------------------------------------

		template<class P, class C, class R = decltype(std::declval<Q>() + std::declval<P>())>
		decltype(auto) operator+(const pmatrix<P, C>& m) const
		{
			if (Dim() == m.Dim())
			{
				TMatrix<R> res(Dim());

				for (unsigned idx = 0; idx < res._data.size(); idx++)
				{
					res._data[idx] = (*this)(idx) + m(idx);
				}

				return res;
			}
			else
				throw IncompatibleDimensions();
		}

		template<class P, class C>
		void operator+=(const pmatrix<P, C>& m)
		{
			if (Dim() == m.Dim())
			{
				for (unsigned idx = 0; idx < _container._data.size(); idx++)
				{
					(*this)(idx) += m.at(idx);
				}
			}
			else
				throw IncompatibleDimensions();
		}

		template<class P, class C, class R = decltype(std::declval<Q>() - std::declval<P>())>
		decltype(auto) operator-(const pmatrix<P, C>& m) const
		{
			if (Dim() == m.Dim())
			{
				TMatrix<R> res(Dim());

				for (unsigned idx = 0; idx < res._data.size(); idx++)
				{
					res._data[idx] = this->at(idx) - m.at(idx);
				}

				return res;
			}
			else
				throw IncompatibleDimensions();
		}

		template<class P, class C>
		void operator-=(const pmatrix<P, C>& m)
		{
			if (Dim() == m.Dim())
			{
				for (unsigned idx = 0; idx < _container._data.size(); idx++)
				{
					(*this)(idx) -= m.at(idx);
				}
			}
			else
				throw IncompatibleDimensions();
		}

		decltype(auto) operator*(const Q& a) const
		{
			TMatrix res(Dim());

			for (unsigned idx = 0; idx < res._data.size(); idx++)
			{
				res._data[idx] = this->at(idx) * a;
			}

			return res;
		}

		void operator*=(const Q& a)
		{
			for (auto& e : _container._data)
			{
				e *= a;
			}
		}

		template<class P, class C, class R = decltype(std::declval<Q>() * std::declval<P>())>
		decltype(auto) operator*(const pmatrix<P, C>& m) const
		{
			if (NCols() == m.NRows())
			{
				TMatrix<R> res(NRows(), m.NCols());

				unsigned idx = 0;
				for (const auto& row : rows())
				{
					for (const auto& col: m.cols())
					{
						res._data[idx++] = row & col;
					}
				}

				return res;
			}
			else
				throw IncompatibleDimensions();
		}

		/*template<class P, class C>
		void operator*=(const pmatrix<P, C>& m)
		{
		}*/

		decltype(auto) operator/(const Q& a) const
		{
			TMatrix<Q> res(Dim());

			for (unsigned idx = 0; idx < res._data.size(); idx++)
			{
				res._data[idx] = this->at(idx) / a;
			}

			return res;
		}

		void operator/=(const Q& a)
		{
			for (auto& e : _container._data)
			{
				e /= a;
			}
		}

		TMatrix<Q> operator!() const       ///invert
		{
			const double err = 1e-15;

			if(NRows() == NCols())
			{
				unsigned N = NRows();
				T _t;

				TMatrix<Q> A = *this;
				TMatrix<Q> I(Dim(), [](auto row, auto col){ return row == col; });

				for(unsigned i = 0; i < N; i++)
				{
					const auto& a0 = A[i][i];

					if(std::abs(a0) < err)
					{
						for(unsigned b = i+1; b < N; b++)
						{
							if(std::abs(A[b][i]) > err)
							{
								A[i] += A[b];
								I[i] += I[b];

								/*for(unsigned t = 0; t < N; t++)
								{
									A[i][t] += A[b][t];
									I[i][t] += I[b][t];
								}*/
								break;
							}
						}
					}

					if(std::abs(a0) < err)
						throw SingularMatrix();
						///throw "Matrix doesn\'t have an inverse one";

					if(a0 != 1.0)
					{
						_t = a0;

						for(unsigned x = 0; x < N; x++)
						{
							A[i][x] = A[i][x] /_t;
							I[i][x] = I[i][x] /_t;
						}
					}

					for(unsigned y = 0; y < N; y++)
					{
						_t = A[y][i];

						for(unsigned x = 0; x < N; x++)
						{
							if(y != i && _t != 0.0)
							{
								A[y][x] -= A[i][x] *_t;
								I[y][x] -= I[i][x] *_t;
							}
						}
					}
				}

				return I;
			}
			else
				throw NonSquareMatrix();
		}
	};

    ///--------------------------------------------------------

	std::vector<T> _data;
	unsigned _rows, _cols;

    public:

    TMatrix(): _rows(0), _cols(0)
    {
    }

    TMatrix(unsigned rows, unsigned cols): _data(rows*cols), _rows(rows), _cols(cols)
    {
    }

	TMatrix(unsigned rows, unsigned cols, const T& val) : _data(rows*cols), _rows(rows), _cols(cols)
	{
		for (unsigned idx = 0; idx < _data.size(); idx++)
		{
			_data[idx] = val;
		}
	}

	TMatrix(const std::pair<unsigned, unsigned>& dim) : _data(std::get<0>(dim)*std::get<1>(dim)), _rows(std::get<0>(dim)), _cols(std::get<1>(dim))
	{
	}

	TMatrix(const std::pair<unsigned, unsigned>& dim, const T& val) : _data(std::get<0>(dim)*std::get<1>(dim)), _rows(std::get<0>(dim)), _cols(std::get<1>(dim))
	{
		for (unsigned idx = 0; idx < _data.size(); idx++)
		{
			_data[idx] = val;
		}
	}

	template<class TFunction, class = decltype(std::declval<TFunction>()(0, 0))>
	TMatrix(unsigned rows, unsigned cols, const TFunction& f) : _data(rows*cols), _rows(rows), _cols(cols)
	{
		for (unsigned row = 0; row < _rows; row++)
		{
			for (unsigned col = 0; col < _cols; col++)
			{
				_data[row*_cols + col] = f(row, col);
			}
		}
	}

	template<class TFunction, class = decltype(std::declval<TFunction>()(0, 0))>
	TMatrix(const std::pair<unsigned, unsigned>& dim, const TFunction& f) : _data(std::get<0>(dim)*std::get<1>(dim)), _rows(std::get<0>(dim)), _cols(std::get<1>(dim))
	{
		for (unsigned row = 0; row < _rows; row++)
		{
			for (unsigned col = 0; col < _cols; col++)
			{
				_data[row*_cols + col] = f(row, col);
			}
		}
	}

    TMatrix(const std::vector<std::vector<T>>& data): _rows(0), _cols(0)
    {
        if (data.size())
        {
            _rows = data.size();
            _cols = data.begin()->size();

            _data.resize(_rows*_cols);

            unsigned idx = 0;
            for (const auto& row: data)
            {
                for (const auto& e: row)
                {
                    _data[idx++] = e;
                }
            }
        }
    }

    TMatrix(const std::initializer_list<std::vector<T>>& data): _rows(0), _cols(0)
    {
        if (data.size())
        {
            _rows = data.size();
            _cols = data.begin()->size();

            _data.resize(_rows*_cols);

            unsigned idx = 0;
            for (const auto& row: data)
            {
                for (const auto& e: row)
                {
                    _data[idx++] = e;
                }
            }
        }
    }

    ~TMatrix()
    {
    }

    unsigned NRows() const
    {
        return _rows;
    }

    unsigned NCols() const
    {
        return _cols;
    }

    void resize(unsigned rows, unsigned cols)
    {
        _rows = rows;
        _cols = cols;

        _data.resize(_rows*_cols);
    }

	void resize(const std::pair<unsigned, unsigned>& dim)
	{
		resize(std::get<0>(dim), std::get<1>(dim));
	}

	decltype(auto) rows()
	{
		return pmatrix<T, TMatrix<T>>(*this, _data.data(), _data.data() + _data.size(), NCols(), 1, false);
	}

	decltype(auto) cols()
	{
		return pmatrix<T, TMatrix<T>>(*this, _data.data(), _data.data() + NCols(), 1, NCols(), true);
	}

	decltype(auto) rows() const
	{
		return pmatrix<T, const TMatrix<T>>(*this, _data.data(), _data.data() + _data.size(), NCols(), 1, false);
	}

	decltype(auto) cols() const
	{
		return pmatrix<T, const TMatrix<T>>(*this, _data.data(), _data.data() + NCols(), 1, NCols(), true);
	}

	operator pmatrix<T, const TMatrix<T>>() const
	{
		return rows();
	}

    pvector<T> operator[](int row)
    {
        return rows()[row];
    }

    ///-----------------------------------------------------

	std::pair<unsigned, unsigned> Dim() const
	{
		return {_rows, _cols};
	}

	explicit operator const T&() const
	{
		if (_data.size() == 1)
			return _data[0];
		else
			throw InvalidScalarCast();
	}

	template<class MType>
	void operator+=(const MType& m)
	{
		rows() += m.rows();
	}

	template<class MType>
	void operator-=(const MType& m)
	{
		rows() -= m.rows();
	}

	decltype(auto) operator*(const T& a) const
	{
		return rows()*a;
	}

	void operator*=(const T& a)
	{
		rows() *= a;
	}

	decltype(auto) operator/(const T& a) const
	{
		return rows()/a;
	}

	void operator/=(const T& a)
	{
		rows() /= a;
	}

	T Sqr() const
    {
        T res = 0;

        for (const auto& e: _data)
		{
			res += e*e;
		}

        return res;
    }

    T Norm() const
    {
        return _data.size() ? sqrt(Sqr()/_data.size()) : 0;
    }

	TMatrix<T> operator!() const       ///invert
    {
        return !rows();
    }

	/// scalar product of matrices treated as one-dimensional arrays
	T operator,(const TMatrix<T>& m) const
	{
		if (_data.size() == m._data.size())
		{
			T res = 0;

			for (unsigned idx = 0; idx < _data.size(); idx++)
			{
				res += _data[idx]*m._data[idx];
			}

			return res;
		}
		else
			throw IncompatibleDimensions();
	}

    /*
    TMatrix operator*() const    ///transpose
    {
        TMatrix v;

        v.Set(M, N);

        for(unsigned x=0;x<N;x++)
        {
            for(unsigned y=0;y<M;y++)
            {
                v.a[y][x] = a[x][y];
            }
        }

        return v;
    }

    ///-----------------------------------------------------

    TMatrix Cov() const
    {
    	TMatrix v = *this;

    	if(this->N == 1) v = v*(*v);
    	else if(this->M == 1) v = (*v)*(v);

    	return v;
    }

    T Trace() const
    {
    	T res = 0;

    	if(a)
    	{
    		for(unsigned x=0;x<_Min(this->N, this->M);x++)
    		{
				res += a[x][x];
    		}
    	}

    	return res;
    }

    TMatrix Tr() const
    {
    	TMatrix v = TMatrix(1, _Min(this->N, this->M));

    	if(a)
    	{
    		for(unsigned x=0;x<_Min(this->N, this->M);x++)
    		{
				v.a[0][x] = a[x][x];
    		}
    	}

    	return v;
    }

	T Sum() const
	{
		T res = 0;

		if(a)
		{
			for(unsigned n=0;n<this->N;n++)
			{
				for(unsigned m=0;m<this->M;m++)
				{
					res += a[n][m];
				}
			}
		}

		return res;
	}

    ///-----------------------------------------------------

	template<class Q>
	friend TMatrix<T> operator^(Q&& f, const TMatrix<T>& m)
	{
		TMatrix<T> v(m);

		if(m.a)
		{
			for(unsigned x=0; x<m.N; x++)
			{
				for(unsigned y=0; y<m.M; y++)
				{
					v.a[x][y] = f(m.a[x][y]);
				}
			}
		}

		return v;
	}

    template<class Q, class... Args>
    TMatrix<T>& apply(Q&& f, Args&&... args)
    {
        if(a)
        {
            for(unsigned x=0; x<N; x++)
            {
                for(unsigned y=0; y<M; y++)
                {
                    f(a[x][y], args...);
                }
            }
        }

        return *this;
    }

    TMatrix<T>& apply(T (*f)(T))
    {
        if(a)
        {
            for(unsigned x=0; x<N; x++)
            {
                for(unsigned y=0; y<M; y++)
                {
                    a[x][y] = f(a[x][y]);
                }
            }
        }

        return *this;
    }*/

	template<class Q>
	friend Q mean(const TMatrix<Q>& m)
	{
		Q res = 0;

		for (const auto& e : m._data)
		{
			res += e;
		}

		return m._data.size() ? res/m._data.size() : 0;
	}

	std::string ToStr(const StrFormat& format = Plain) const
	{
		std::string res;
		char s[32];

		switch(format)
		{
			case List:

				for (const auto& row : rows())
				{
					if (!res.empty())
						res += ", ";

					res += '{';

					bool _tab = false;
					for (const auto& e : row)
					{
						if (_tab)
						{
							res += ", ";
						}
						else
							_tab = true;

						sprintf(s, "%g", (double)e);
						res += s;
					}

					res += '}';
				}

				return "{" + res + "}";

			default:	/// plain
				for (const auto& row : rows())
				{
					if (!res.empty())
						res += '\n';

					bool _tab = false;
					for (const auto& e : row)
					{
						if (_tab)
						{
							res += '\t';
						}
						else
							_tab = true;

						res += std::to_string(e);
					}
				}
				break;
		}

		return res;
	}
};

template<class M1, class M2, class = decltype(std::declval<M1>().rows()), class = decltype(std::declval<M2>().rows())>
decltype(auto) operator*(const M1& l, const M2& r)
{
    return l.rows() * r.rows();
}

template<class M1, class M2, class = decltype(std::declval<M1>().rows()), class = decltype(std::declval<M2>().rows())>
decltype(auto) operator+(const M1& l, const M2& r)
{
    return l.rows() + r.rows();
}

template<class M1, class M2, class = decltype(std::declval<M1>().rows()), class = decltype(std::declval<M2>().rows())>
decltype(auto) operator-(const M1& l, const M2& r)
{
    return l.rows() - r.rows();
}

///-----------------------------------------------------------------------------

typedef TMatrix<double> matrix;

///-----------------------------------------------------------------------------

template<class PMatrix>
decltype(auto) operator*(const typename PMatrix::type& a, const PMatrix& m)
{
	TMatrix<typename PMatrix::type> res(m.Dim());

	for (unsigned idx = 0; idx < res._data.size(); idx++)
	{
		res._data[idx] = a * m.rows().at(idx);
	}

	return res;
}

#if defined(_GLIBCXX_IOSTREAM) || defined(_IOSTREAM_)

//#define COLOURED_NUMBERS
#define PRETTY_TABLE

template<class PMatrix>
inline auto operator<<(std::ostream& out, PMatrix&& m) -> decltype(m.rows(), out)
{
	#ifdef COLOURED_NUMBERS

	const auto WHITE = FOREGROUND_BLUE | FOREGROUND_GREEN | FOREGROUND_RED;

	#endif

	#ifdef PRETTY_TABLE

	out.flags(std::ios::left | std::ios::right | std::ios::showpos);
	std::cout.setf(std::ios::fixed);//, std::ios::floatfield);

	#endif

	for (const auto& row : m.rows())
	{
		bool _tab = false;
		for (const auto& e : row)
		{
			if (_tab)
			{
				#ifdef PRETTY_TABLE

				out << '\t';

				#else

				out << ' ';

				#endif
			}
			else
				_tab = true;

			#ifdef COLOURED_NUMBERS

			if (e < 0)
				SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), FOREGROUND_BLUE | FOREGROUND_INTENSITY);

			#endif

			out << e;

			#ifdef COLOURED_NUMBERS

			SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), WHITE);

			#endif
		}
		out << '\n';
	}

	#ifdef PRETTY_TABLE

	std::cout.unsetf(std::ios::fixed);
	std::cout.unsetf(std::ios::left | std::ios::right | std::ios::showpos);

	#endif

	return out;
}

#endif

///-------------------------------------------------------------------

/// Linear systems solver (conjugate gradients method)
///
///ret x: |A*x - b| < C*dr
///
template<class MType>
inline MType GrLSolve(const MType& A, const MType& b, const double dr = 1e-12)
{
    if(A.NRows() == b.NRows())
    {
		unsigned n = 0;

		MType x(A.NCols(), 1, 0);
		MType p(A.NCols(), 1, 0);

		double R0 = -1, R = 0;

        MType q, r, dx;
		auto tA = A.cols();

        while(n <= 2*A.NCols())
        {
            if(n) r -= q/(p, q);
            else
            {
                r = tA*(A*x - b);

                R0 = r.Sqr();

                if(R0 == 0) break;
            }

            R = r.Sqr();

            p += r/R;
			q = tA*(A*p);

            dx = p/(p, q);
			x -= dx;

            n++;

            if(R <= dr*dr*R0) break;
        }

        //printf("\n[%d/%d iterations: \263*A\371(A\371x - b)\263 = %g]\n\n", n, 2*A.NCols(), sqrt(R));

		return x;
    }
	else
		throw IncompatibleDimensions();
}

#endif // TMATRIX_H_INCLUDED
