// Formula.h
// Copyright (C) 2011 Miguel Vargas (miguel.vargas@gmail.com)
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

#ifndef _Formula_h_
#define _Formula_h_

#include <Basic/Assert.h>
#include <Basic/Debug.h>
#include <Basic/Memory.h>
#include <Basic/String.h>
#include <Container/List.h>

#include <math.h>
#include <stdarg.h>
#include <string.h>

// Examples:
//  exp(-(t - 5)*(t - 5)/0.1); t
//  sqrt(x*x + y*y); x; y
//
// Functions:
//  abs    acos  asin   atan
//  atan2  ceil  cos    cosh
//  tan    exp   floor  log
//  log10  mod   pow    sin
//  sinh   sqrt  tan    tanh


template <typename T>
class Formula
{
	public:

		class Exception : public ::Exception {};

		class ExceptionBracket : public Exception {};

		class ExceptionSyntax : public Exception {};

		class ExceptionUnknownFunction : public Exception {};

		class ExceptionVariableName : public Exception {};


		static Exception exception;

		static ExceptionBracket exception_bracket;

		static ExceptionSyntax exception_syntax;

		static ExceptionUnknownFunction exception_unknown_function;

		static ExceptionVariableName exception_variable_name;


		Formula() throw()
		:	root_item((FormulaItem*)0),
			items(),
			variables()
		{
		}


		Formula(const char* text) throw(Memory::Exception, ExceptionBracket, ExceptionSyntax, ExceptionUnknownFunction, ExceptionVariableName)
		:	root_item((FormulaItem*)0),
			items(),
			variables()
		{
			Parse(text);
		}


		void Define(const char* text) throw(Memory::Exception, ExceptionBracket, ExceptionSyntax, ExceptionUnknownFunction, ExceptionVariableName)
		{
			Assert(text);

			try
			{
				Clean();

				String formula_copy((int)strlen(text) + 1);
				const char* pre_formula = text;

				// Count and split variables
				int variables_count = 0;
				char* variable_names = (char*)0;
				for (char* f = formula_copy.data; ; ++pre_formula)
				{
					if (*pre_formula == ';')
					{
						*f = '\0';
						++f;
						if (variables_count == 0)
						{
							variable_names = f;
						}
						++variables_count;
					}
					else if (*pre_formula == '\0')
					{
						*f = '\0';
						break;
					}
					else
					{
						*f = *pre_formula;
						++f;
					}
				}

				// Create variables
				for (int v = 1; v <= variables_count; ++v, ++variable_names)
				{
					char* name = ParseVariableDefinition(variable_names);
					Variable* new_variable = new Variable(name);
					if (!new_variable)
					{
						Throw(Memory::exception);
					}
					variables.AppendLast(new_variable);
				}

				char* formula = formula_copy.data;
				root_item = Parse(formula);
			}
			catch (Exception&)
			{
				ReThrow();
			}
		}


		~Formula() throw()
		{
			for (ListItem<FormulaItem*>* list_item = items.first; list_item; list_item = list_item->next)
			{
				delete list_item->value;
			}
			for (ListItem<Variable*>* list_item = variables.first; list_item; list_item = list_item->next)
			{
				delete list_item->value;
			}
		}


		void Clean() throw()
		{
			for (ListItem<FormulaItem*>* list_item = items.first; list_item; list_item = list_item->next)
			{
				delete list_item->value;
			}
			for (ListItem<Variable*>* list_item = variables.first; list_item; list_item = list_item->next)
			{
				delete list_item->value;
			}
			root_item = (FormulaItem*)0;
		}


		T operator () (void) throw()
		{
			Assert(root_item);
			return root_item->Evaluate();
		}


		T operator () (T v1, ...) throw()
		{
			Assert(root_item);

			if (variables.size > 0)
			{
				variables.first->value->var = v1;

				va_list variadic;
				va_start(variadic, v1);
				for (register ListItem<Variable*>* item = variables.first->next; item; item = item->next)
				{
					item->value->var = va_arg(variadic, T);
				}
				va_end(variadic);
			}
			return root_item->Evaluate();
		}


	private:

		struct FormulaItem
		{
			virtual ~FormulaItem() throw()
			{
			}

			virtual T Evaluate(void) const throw() = 0;
		};


		struct Variable : FormulaItem
		{
			String name;
			T var;

			Variable(const char* name) throw()
			:	FormulaItem(),
				name(name),
				var()
			{
				Assert(name);
			}

			Variable& operator = (const Variable& variable) throw()
			{
				name = variable.name;
				var = variable.var;
				return *this;
			}

			virtual T Evaluate(void) const throw()
			{
				return var;
			}
		};


		struct Constant : FormulaItem
		{
			T value;

			Constant(T value) throw()
			:	FormulaItem(),
				value(value)
			{
			}

			virtual T Evaluate(void) const throw()
			{
				return value;
			}
		};


		struct UnaryItem : FormulaItem
		{
			FormulaItem* a;

			UnaryItem(FormulaItem* a) throw()
			:	FormulaItem(),
				a(a)
			{
				Assert(a);
			}
		};


		struct BinaryItem : FormulaItem
		{
			FormulaItem* a;
			FormulaItem* b;

			BinaryItem(FormulaItem* a, FormulaItem* b) throw()
			:	FormulaItem(),
				a(a),
				b(b)
			{
				Assert(a);
				Assert(b);
			}
		};


		struct Minus : UnaryItem
		{
			Minus(FormulaItem* a) throw()
			:	UnaryItem(a)
			{
			}

			virtual T Evaluate(void) const throw()
			{
				return -this->a->Evaluate();
			}
		};


		struct Addition : BinaryItem
		{
			Addition(FormulaItem* a, FormulaItem* b) throw()
			:	BinaryItem(a, b)
			{
			}

			virtual T Evaluate(void) const throw()
			{
				return this->a->Evaluate() + this->b->Evaluate();
			}
		};


		struct Subtraction : BinaryItem
		{
			Subtraction(FormulaItem* a, FormulaItem* b) throw()
			:	BinaryItem(a, b)
			{
			}

			virtual T Evaluate(void) const throw()
			{
				return this->a->Evaluate() - this->b->Evaluate();
			}
		};


		struct Multiplication : BinaryItem
		{
			Multiplication(FormulaItem* a, FormulaItem* b) throw()
			:	BinaryItem(a, b)
			{
			}

			virtual T Evaluate(void) const throw()
			{
				return this->a->Evaluate()*this->b->Evaluate();
			}
		};


		struct Division : BinaryItem
		{
			Division(FormulaItem* a, FormulaItem* b) throw()
			:	BinaryItem(a, b)
			{
			}

			virtual T Evaluate(void) const throw()
			{
				return this->a->Evaluate()/this->b->Evaluate();
			}
		};


		struct Absolute : UnaryItem
		{
			Absolute(FormulaItem* a) throw()
			:	UnaryItem(a)
			{
			}

			virtual T Evaluate(void) const throw()
			{
				return fabs(this->a->Evaluate());
			}
		};


		struct ArcCosine : UnaryItem
		{
			ArcCosine(FormulaItem* a) throw()
			:	UnaryItem(a)
			{
			}

			virtual T Evaluate(void) const throw()
			{
				return acos(this->a->Evaluate());
			}
		};


		struct ArcSine : UnaryItem
		{
			ArcSine(FormulaItem* a) throw()
			:	UnaryItem(a)
			{
			}

			virtual T Evaluate(void) const throw()
			{
				return asin(this->a->Evaluate());
			}
		};


		struct ArcTangent : UnaryItem
		{
			ArcTangent(FormulaItem* a) throw()
			:	UnaryItem(a)
			{
			}

			virtual T Evaluate(void) const throw()
			{
				return atan(this->a->Evaluate());
			}
		};


		struct ArcTangent2 : BinaryItem
		{
			ArcTangent2(FormulaItem* a, FormulaItem* b) throw()
			:	BinaryItem(a, b)
			{
			}

			virtual T Evaluate(void) const throw()
			{
				return atan2(this->a->Evaluate(), this->b->Evaluate());
			}
		};


		struct Ceil : UnaryItem
		{
			Ceil(FormulaItem* a) throw()
			:	UnaryItem(a)
			{
			}

			virtual T Evaluate(void) const throw()
			{
				return ceil(this->a->Evaluate());
			}
		};


		struct Cosine : UnaryItem
		{
			Cosine(FormulaItem* a) throw()
			:	UnaryItem(a)
			{
			}

			virtual T Evaluate(void) const throw()
			{
				return cos(this->a->Evaluate());
			}
		};


		struct CosineHyperbolic : UnaryItem
		{
			CosineHyperbolic(FormulaItem* a) throw()
			:	UnaryItem(a)
			{
			}

			virtual T Evaluate(void) const throw()
			{
				return cosh(this->a->Evaluate());
			}
		};


		struct Exponential : UnaryItem
		{
			Exponential(FormulaItem* a) throw()
			:	UnaryItem(a)
			{
			}

			virtual T Evaluate(void) const throw()
			{
				return exp(this->a->Evaluate());
			}
		};


		struct Floor : UnaryItem
		{
			Floor(FormulaItem* a) throw()
			:	UnaryItem(a)
			{
			}

			virtual T Evaluate(void) const throw()
			{
				return floor(this->a->Evaluate());
			}
		};


		struct Logarithm : UnaryItem
		{
			Logarithm(FormulaItem* a) throw()
			:	UnaryItem(a)
			{
			}

			virtual T Evaluate(void) const throw()
			{
				return log(this->a->Evaluate());
			}
		};


		struct Logarithm10 : UnaryItem
		{
			Logarithm10(FormulaItem* a) throw()
			:	UnaryItem(a)
			{
			}

			virtual T Evaluate(void) const throw()
			{
				return log10(this->a->Evaluate());
			}
		};


		struct Modulo : BinaryItem
		{
			Modulo(FormulaItem* a, FormulaItem* b) throw()
			:	BinaryItem(a, b)
			{
			}

			virtual T Evaluate(void) const throw()
			{
				return fmod(this->a->Evaluate(), this->b->Evaluate());
			}
		};


		struct Power : BinaryItem
		{
			Power(FormulaItem* a, FormulaItem* b) throw()
			:	BinaryItem(a, b)
			{
			}

			virtual T Evaluate(void) const throw()
			{
				return pow(this->a->Evaluate(), this->b->Evaluate());
			}
		};


		struct Sine : UnaryItem
		{
			Sine(FormulaItem* a) throw()
			:	UnaryItem(a)
			{
			}

			virtual T Evaluate(void) const throw()
			{
				return sin(this->a->Evaluate());
			}
		};


		struct SineHyperbolic : UnaryItem
		{
			SineHyperbolic(FormulaItem* a) throw()
			:	UnaryItem(a)
			{
			}

			virtual T Evaluate(void) const throw()
			{
				return sinh(this->a->Evaluate());
			}
		};


		struct SquareRoot : UnaryItem
		{
			SquareRoot(FormulaItem* a) throw()
			:	UnaryItem(a)
			{
			}

			virtual T Evaluate(void) const throw()
			{
				return sqrt(this->a->Evaluate());
			}
		};


		struct Tangent : UnaryItem
		{
			Tangent(FormulaItem* a) throw()
			:	UnaryItem(a)
			{
			}

			virtual T Evaluate(void) const throw()
			{
				return tan(this->a->Evaluate());
			}
		};


		struct TangentHyperbolic : UnaryItem
		{
			TangentHyperbolic(FormulaItem* a) throw()
			:	UnaryItem(a)
			{
			}

			virtual T Evaluate(void) const throw()
			{
				return tanh(this->a->Evaluate());
			}
		};


		FormulaItem* Parse(char*& formula) throw(Memory::Exception, ExceptionBracket, ExceptionSyntax, ExceptionUnknownFunction)
		{
			FormulaItem* base_item = ParseMonomial(formula);

			for ( ; ; )
			{
				while ((*formula == ' ') || (*formula == '\t'))
				{
					++formula;
				}

				FormulaItem* new_item = (FormulaItem*)0;
				if (*formula == '+')
				{
					++formula;
					new_item = new Addition(base_item, ParseMonomial(formula));
				}
				else if (*formula == '-')
				{
					++formula;
					new_item = new Subtraction(base_item, ParseMonomial(formula));
				}
				else if (*formula == '\0')
				{
					break;
				}
				else
				{
					Throw(Formula::exception_syntax);
				}

				if (!new_item)
				{
					Throw(Memory::exception);
				}
				items.AppendLast(new_item);
				base_item = new_item;
			}
			return base_item;
		}


		FormulaItem* ParseMonomial(char*& formula) throw(Memory::Exception, ExceptionBracket, ExceptionSyntax, ExceptionUnknownFunction)
		{
			FormulaItem* base_item = ParseItem(formula);

			for ( ; ; )
			{
				while ((*formula == ' ') || (*formula == '\t'))
				{
					++formula;
				}

				FormulaItem* new_item = (FormulaItem*)0;
				if (*formula == '*')
				{
					++formula;
					new_item = new Multiplication(base_item, ParseItem(formula));
				}
				else if (*formula == '/')
				{
					++formula;
					new_item = new Division(base_item, ParseItem(formula));
				}
				else if ((*formula == '+') || (*formula == '-') || (*formula == '\0'))
				{
					break;
				}
				else
				{
					Throw(Formula::exception_syntax);
				}

				if (!new_item)
				{
					Throw(Memory::exception);
				}
				items.AppendLast(new_item);
				base_item = new_item;
			}
			return base_item;
		}


		FormulaItem* ParseItem(char*& formula) throw(Memory::Exception, ExceptionBracket, ExceptionSyntax, ExceptionUnknownFunction)
		{
			while ((*formula == ' ') || (*formula == '\t'))
			{
				++formula;
			}

			if ((*formula >= '0') && (*formula <= '9'))
			{
				return ParseNumber(formula);
			}
			else if (((*formula >= 'a') && (*formula <= 'z')) || ((*formula >= 'A') && (*formula <= 'Z')))
			{
				return ParseSymbol(formula);
			}
			else if (*formula == '(')
			{
				return ParseBrackets(formula);
			}
			else if (*formula == '+')
			{
				++formula;
				return ParseItem(formula);
			}
			else if (*formula == '-')
			{
				++formula;
				FormulaItem* new_item = new Minus(ParseItem(formula));
				if (!new_item)
				{
					Throw(Memory::exception);
				}
				items.AppendLast(new_item);
				return new_item;
			}

			Throw(Formula::exception_syntax);
		}


		FormulaItem* ParseBrackets(char*& formula) throw(Memory::Exception, ExceptionBracket, ExceptionSyntax, ExceptionUnknownFunction)
		{
			++formula;
			char* sub_formula = formula;
			int bracket_level = 1;
			for ( ; *formula; ++formula)
			{
				if (*formula == ')')
				{
					--bracket_level;
					if (bracket_level == 0)
					{
						*formula = '\0';
						++formula;
						return Parse(sub_formula);
					}
				}
				else if (*formula == '(')
				{
					++bracket_level;
				}
			}
			Throw(Formula::exception_bracket;);
		}


		FormulaItem* ParseBracketsWithComma(char*& formula, FormulaItem*& first_item) throw(Memory::Exception, ExceptionBracket, ExceptionSyntax, ExceptionUnknownFunction)
		{
			first_item = (FormulaItem*)0;

			++formula;
			char* sub_formula = formula;
			int bracket_level = 1;
			for ( ; *formula; ++formula)
			{
				if (*formula == ',')
				{
					if (bracket_level == 1)
					{
						if (first_item)
						{
							Throw(Formula::exception_syntax);
						}

						*formula = '\0';
						++formula;
						first_item = Parse(sub_formula);
						sub_formula = formula;
					}
				}
				if (*formula == ')')
				{
					--bracket_level;
					if (bracket_level == 0)
					{
						if (!first_item)
						{
							Throw(Formula::exception_syntax);
						}

						*formula = '\0';
						++formula;
						return Parse(sub_formula);
					}
				}
				else if (*formula == '(')
				{
					++bracket_level;
				}
			}
			Throw(Formula::exception_bracket;);
		}


		FormulaItem* ParseSymbol(char*& formula) throw(Memory::Exception, ExceptionBracket, ExceptionSyntax, ExceptionUnknownFunction)
		{
			char* symbol = formula;
			int length = 0;
			do
			{
				++length;
				++formula;
			} while (((*formula >= 'a') && (*formula <= 'z')) || ((*formula >= 'A') && (*formula <= 'Z')) || ((*formula >= '0') && (*formula <= '9')) || (*formula == '_'));

			for (ListItem<Variable*>* item = variables.first; item; item = item->next)
			{
				if (strncmp(symbol, item->value->name.data, length) == 0)
				{
					return item->value;
				}
			}

			while ((*formula == ' ') || (*formula =='\t'))
			{
				*formula = '\0';
				++formula;
			}

			if (*formula == '(')
			{
				FormulaItem* first_item = (FormulaItem*)0;
				FormulaItem* second_item = (FormulaItem*)0;

				*formula = '\0';
				FormulaItem* new_item = (FormulaItem*)0;
				if (strcmp(symbol, "abs") == 0)
				{
					new_item = new Absolute(ParseBrackets(formula));
				}
				else if (strcmp(symbol, "acos") == 0)
				{
					new_item = new ArcCosine(ParseBrackets(formula));
				}
				else if (strcmp(symbol, "asin") == 0)
				{
					new_item = new ArcSine(ParseBrackets(formula));
				}
				else if (strcmp(symbol, "atan") == 0)
				{
					new_item = new ArcTangent(ParseBrackets(formula));
				}
				else if (strcmp(symbol, "atan2") == 0)
				{
					second_item = ParseBracketsWithComma(formula, first_item);
					new_item = new ArcTangent2(first_item, second_item);
				}
				else if (strcmp(symbol, "ceil") == 0)
				{
					new_item = new Ceil(ParseBrackets(formula));
				}
				else if (strcmp(symbol, "cos") == 0)
				{
					new_item = new Cosine(ParseBrackets(formula));
				}
				else if (strcmp(symbol, "cosh") == 0)
				{
					new_item = new CosineHyperbolic(ParseBrackets(formula));
				}
				else if (strcmp(symbol, "exp") == 0)
				{
					new_item = new Exponential(ParseBrackets(formula));
				}
				else if (strcmp(symbol, "floor") == 0)
				{
					new_item = new Floor(ParseBrackets(formula));
				}
				else if (strcmp(symbol, "log") == 0)
				{
					new_item = new Logarithm(ParseBrackets(formula));
				}
				else if (strcmp(symbol, "log10") == 0)
				{
					new_item = new Logarithm10(ParseBrackets(formula));
				}
				else if (strcmp(symbol, "mod") == 0)
				{
					second_item = ParseBracketsWithComma(formula, first_item);
					new_item = new Modulo(first_item, second_item);
				}
				else if (strcmp(symbol, "pow") == 0)
				{
					second_item = ParseBracketsWithComma(formula, first_item);
					new_item = new Power(first_item, second_item);
				}
				else if (strcmp(symbol, "sin") == 0)
				{
					new_item = new Sine(ParseBrackets(formula));
				}
				else if (strcmp(symbol, "sinh") == 0)
				{
					new_item = new SineHyperbolic(ParseBrackets(formula));
				}
				else if (strcmp(symbol, "sqrt") == 0)
				{
					new_item = new SquareRoot(ParseBrackets(formula));
				}
				else if (strcmp(symbol, "tan") == 0)
				{
					new_item = new Tangent(ParseBrackets(formula));
				}
				else if (strcmp(symbol, "tanh") == 0)
				{
					new_item = new TangentHyperbolic(ParseBrackets(formula));
				}
				else
				{
					Throw(Formula::exception_unknown_function);
				}
				if (!new_item)
				{
					Throw(Memory::exception);
				}
				items.AppendLast(new_item);
				return new_item;
			}
			
			Throw(Formula::exception_syntax);
		}


		FormulaItem* ParseNumber(char*& formula) throw(Memory::Exception, ExceptionSyntax)
		{
			enum
			{
				ns_integer,
				ns_radix_point,
				ns_fractional,
				ns_exp_symbol,
				ns_exp_sign,
				ns_exponent,
				ns_end
			} state = ns_integer;
			static const bool from_to_state[ns_end][ns_end + 1] =
			{
				{true , true , false, true , false, false, true },
				{false, false, true , false, false, false, false},
				{false, false, true , true , false, false, true },
				{false, false, false, false, true , true , false},
				{false, false, false, false, false, true , false},
				{false, false, false, false, false, true , true }
			};

			T integer = 0;
			T fractional = 0;
			T divisor = 1;
			T exp_sign = 1;
			T exponent = 0;
			for ( ; ; ++formula)
			{
				if ((*formula >= '0') && (*formula <= '9'))
				{
					if (from_to_state[state][ns_integer])
					{
						state = ns_integer;
						integer *= 10;
						integer += *formula - '0';
						continue;
					}
					else if (from_to_state[state][ns_fractional])
					{
						state = ns_fractional;
						fractional *= 10;
						fractional += *formula - '0';
						divisor *= 10;
						continue;
					}
					else if (from_to_state[state][ns_exponent])
					{
						state = ns_exponent;
						exponent *= 10;
						exponent += *formula - '0';
						continue;
					}
				}
				else if (*formula == '.')
				{
					if (from_to_state[state][ns_radix_point])
					{
						state = ns_radix_point;
						continue;
					}
				}
				else if ((*formula == 'e') || (*formula == 'E'))
				{
					if (from_to_state[state][ns_exp_symbol])
					{
						state = ns_exp_symbol;
						continue;
					}
				}
				else if (*formula == '+')
				{
					if (from_to_state[state][ns_exp_sign])
					{
						state = ns_exp_sign;
						continue;
					}
					else if (from_to_state[state][ns_end])
					{
						break;
					}

				}
				else if (*formula == '-')
				{
					if (from_to_state[state][ns_exp_sign])
					{
						state = ns_exp_sign;
						exp_sign = -1;
						continue;
					}
					else if (from_to_state[state][ns_end])
					{
						break;
					}
				}
				else
				{
					if (from_to_state[state][ns_end])
					{
						break;
					}
				}
				Throw(Formula::exception_syntax);
			}
			T number = (integer*divisor + fractional)*pow(10, exp_sign*exponent)/divisor;

			FormulaItem* new_item = new Constant(number);
			if (!new_item)
			{
				Throw(Memory::exception);
			}
			items.AppendLast(new_item);

			return new_item;
		}


		char* ParseVariableDefinition(char*& formula) throw()
		{
			enum
			{
				vs_start,
				vs_pre_space,
				vs_var_head,
				vs_var_tail,
				vs_post_space,
				vs_end
			} state = vs_start;
			static const bool from_to_state[vs_end][vs_end + 1] =
			{
				{false, true , true,  false, false, false},
				{false, true , true,  false, false, false},
				{false, false, false, true , true , true },
				{false, false, false, true , true , true },
				{false, false, false, false, true , true }
			};

			char* name = (char*)0;
			for ( ; ; ++formula)
			{
				if ((*formula == ' ') || (*formula == '\t'))
				{
					if (from_to_state[state][vs_pre_space])
					{
						state = vs_pre_space;
						continue;
					}
					else if (from_to_state[state][vs_post_space])
					{
						state = vs_post_space;
						*formula = '\0';
						continue;
					}
				}
				else if (((*formula >= '0') && (*formula <= '9')) || (*formula == '_'))
				{
					if (from_to_state[state][vs_var_tail])
					{
						state = vs_var_tail;
						continue;
					}
				}
				else if (((*formula >= 'a') && (*formula <= 'z')) || ((*formula >= 'A') && (*formula <= 'Z')))
				{
					if (from_to_state[state][vs_var_head])
					{
						state = vs_var_head;
						name = formula;
						continue;
					}
					else if (from_to_state[state][vs_var_tail])
					{
						state = vs_var_tail;
						continue;
					}
				}
				else if (*formula == '\0')
				{
					if (from_to_state[state][vs_end])
					{
						break;
					}
				}
				return (char*)0;
			}
			return name;
		}


		FormulaItem* root_item;
		List<FormulaItem*> items;
		List<Variable*> variables;
};

#endif
