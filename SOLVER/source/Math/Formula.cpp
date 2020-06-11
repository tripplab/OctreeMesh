// Formula.cpp
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

#include <Math/Formula.h>


template <>
Formula<float>::Exception Formula<float>::exception;

template <>
Formula<float>::ExceptionBracket Formula<float>::exception_bracket;

template <>
Formula<float>::ExceptionSyntax Formula<float>::exception_syntax;

template <>
Formula<float>::ExceptionUnknownFunction Formula<float>::exception_unknown_function;

template <>
Formula<float>::ExceptionVariableName Formula<float>::exception_variable_name;

template <>
Formula<double>::Exception Formula<double>::exception;

template <>
Formula<double>::ExceptionBracket Formula<double>::exception_bracket;

template <>
Formula<double>::ExceptionSyntax Formula<double>::exception_syntax;

template <>
Formula<double>::ExceptionUnknownFunction Formula<double>::exception_unknown_function;

template <>
Formula<double>::ExceptionVariableName Formula<double>::exception_variable_name;
