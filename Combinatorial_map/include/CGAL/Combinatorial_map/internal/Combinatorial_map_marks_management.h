// Copyright (c) 2010-2011 CNRS and LIRIS' Establishments (France).
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
#ifndef CGAL_COMBINATORIAL_MAP_MARKS_MANAGEMENT_H
#define CGAL_COMBINATORIAL_MAP_MARKS_MANAGEMENT_H 1

#include <bitset>
#include <boost/thread/shared_mutex.hpp>
#include <boost/thread/mutex.hpp>
#include <CGAL/tags.h>
#include <atomic>

namespace CGAL
{
  namespace internal
  {
  template<typename Concurrent_tag>
  struct Size_type
  { typedef std::size_t type; };

  template<>
  struct Size_type<CGAL::Tag_true>
  { typedef std::atomic_size_t type; };
  }
}

  inline
  bool operator==(const std::atomic_size_t& e1, std::size_t e2)
  { return e1.load()==e2; }
  inline
  bool operator==(std::size_t e1, const std::atomic_size_t& e2)
  { return e1==e2.load(); }

  inline
  bool operator!=(const std::atomic_size_t& e1, std::size_t e2)
  { return !(e1==e2); }
  inline
  bool operator!=(std::size_t e1, const std::atomic_size_t& e2)
  { return !(e1==e2); }

  inline
  bool operator>(const std::atomic_size_t& e1, std::size_t e2)
  { return e1.load()>e2; }
  inline
  bool operator>(std::size_t e1, const std::atomic_size_t& e2)
  { return e1>e2.load(); }

  inline
  bool operator<(const std::atomic_size_t& e1, std::size_t e2)
  { return e1.load()<e2; }
  inline
  bool operator<(std::size_t e1, const std::atomic_size_t& e2)
  { return e1<e2.load(); }

  inline
  bool operator>=(const std::atomic_size_t& e1, std::size_t e2)
  { return e1.load()>=e2; }
  inline
  bool operator>=(std::size_t e1, const std::atomic_size_t& e2)
  { return e1>=e2.load(); }

  inline
  bool operator<=(const std::atomic_size_t& e1, std::size_t e2)
  { return e1.load()<=e2; }
  inline
  bool operator<=(std::size_t e1, const std::atomic_size_t& e2)
  { return e1<=e2.load(); }

  namespace CGAL
  {
    namespace internal
    {
    inline
    void set_for_mark(std::atomic_size_t& e1, std::size_t e2)
    { e1.store(e2); }
    inline
    void set_for_mark(std::size_t& e1, const std::atomic_size_t& e2)
    { e1=e2.load(); }
    inline
    void set_for_mark(std::atomic_size_t& e1, const std::atomic_size_t& e2)
    { e1.store(e2.load()); }
    inline
    void set_for_mark(std::size_t& e1, std::size_t e2)
    { e1=e2; }

    inline
    std::size_t get_for_mark(const std::atomic_size_t& e)
    { return e.load(); }
    inline
    std::size_t get_for_mark(std::size_t e)
    { return e; }

    typedef boost::shared_mutex CGAL_CMap_mutex;
    //typedef boost::mutex CGAL_CMap_mutex;

////////////////////////////////////////////////////////////////////////////////
    template<typename Concurrent_tag_, unsigned int NB_MARKS>
    struct Mark_mutexes
    {
      CGAL_CMap_mutex m_global_mutex; // to protect mark array
      CGAL_CMap_mutex m_mask_mutex; // to protect mask bitset
      CGAL_CMap_mutex m_mark_mutexes[NB_MARKS]; // to protect access to each mark separatively

      CGAL_CMap_mutex* global_mutex() const
      { return const_cast<CGAL_CMap_mutex*>(&m_global_mutex); }

      CGAL_CMap_mutex* mask_mutex() const
      { return const_cast<CGAL_CMap_mutex*>(&m_mask_mutex); }

      CGAL_CMap_mutex* mark_mutex(unsigned int m) const
      { return const_cast<CGAL_CMap_mutex*>(&m_mark_mutexes[m]); }
    };

    template<unsigned int NB_MARKS>
    struct Mark_mutexes<CGAL::Tag_false, NB_MARKS>
    {
      CGAL_CMap_mutex* global_mutex() const
      { return nullptr; }

      CGAL_CMap_mutex* mask_mutex() const
      { return nullptr; }

      CGAL_CMap_mutex* mark_mutex(unsigned int /*m*/) const
      { return nullptr; }
    };
////////////////////////////////////////////////////////////////////////////////    

    template<typename Concurrent_tag_>
    struct Mark_array_r_protector: public
        //boost::unique_lock<CGAL_CMap_mutex >
        boost::shared_lock<CGAL_CMap_mutex>
    {
      Mark_array_r_protector(CGAL_CMap_mutex* m):
        //boost::unique_lock<CGAL_CMap_mutex >(*m)
        boost::shared_lock<CGAL_CMap_mutex>(*m)
      {}
    };

    template<>
    struct Mark_array_r_protector<CGAL::Tag_false>
    {
      Mark_array_r_protector(CGAL_CMap_mutex* /*m*/)
      {}
    };
////////////////////////////////////////////////////////////////////////////////    
    template<typename Concurrent_tag_>
    struct Mark_array_w_protector: public boost::unique_lock<CGAL_CMap_mutex >
    {
      Mark_array_w_protector(CGAL_CMap_mutex* m): boost::unique_lock<CGAL_CMap_mutex>(*m)
      {}
    };

    template<>
    struct Mark_array_w_protector<CGAL::Tag_false>
    {
      Mark_array_w_protector(CGAL_CMap_mutex* /*m*/)
      {}
    };
////////////////////////////////////////////////////////////////////////////////    

} //namespace internal

} //namespace CGAL

#endif // CGAL_COMBINATORIAL_MAP_MARKS_MANAGEMENT_H //
// EOF //
