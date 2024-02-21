// Copyright (c) 2024 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
////////////////////////////////////////////////////////////////////////////////
#ifndef LOGGER_H
#define LOGGER_H

#include <iostream>
#include <unistd.h>

class Logger
{
public:
  Logger(const char* f) : mf(f)
  {
    display_msg("BEGIN");
  }

  ~Logger()
  {
    display_msg("END");
  }

  void step(const char* step)
  {
    display_msg(step);
  }

protected:
  void display_msg(const char* m)
  {
    //#ifdef LOG_DISTRIBUTED_LCC
    std::cout<<"{"<<getpid()<<"}"<<"["<<mf<<"] "<<m<<"."<<std::endl;
    //#endif
  }

protected:
  const char* mf;
};

#endif // LOGGER_H
