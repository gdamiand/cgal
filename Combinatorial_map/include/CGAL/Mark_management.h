// Copyright (c) 2026 CNRS and LIRIS' Establishments (France).
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
#ifndef CGAL_MARK_MANAGEMENT_H
#define CGAL_MARK_MANAGEMENT_H 1

#include <CGAL/assertions.h>
#include <bitset>

namespace CGAL {

  template<typename T> // TODO what is T ??
  class Mark_management_bitset_on_dart
  {
  public:
    using size_type=typename T::size_type;
    using Dart_descriptor=typename T::Dart_descriptor;
    using Dart_const_descriptor=typename T::Dart_const_descriptor;

    /// Number of marks
    static const size_type NB_MARKS = 32;

    void reset()
    { mmask_marks.reset(); }

    void swap(Mark_management_bitset_on_dart<T>& other)
    { std::swap(mmask_marks,other.mmask_marks); }

    /// init_dart_mark is called when a new dart is created.
    void init_dart_mark(const T& storage,
                        Dart_descriptor ADart)
    { set_dart_marks(storage, ADart, mmask_marks); }

    /// on_erase_dart is called when a dart is deleted
    void on_erase_dart(const T& storage,
                       Dart_descriptor ADart)
    {}

    /// Set simultaneously all the marks of this dart to a given value.
    void set_dart_marks(const T& storage,
                        Dart_const_descriptor ADart,
                        const std::bitset<NB_MARKS>& amarks) const
    {
      CGAL_assertion( ADart!=nullptr );
      ADart->set_marks(amarks);
    }

    /// Return all the marks of a dart.
    std::bitset<NB_MARKS> get_dart_marks(const T& storage,
                        Dart_const_descriptor ADart) const
    {
      CGAL_assertion( ADart!=nullptr );
      return ADart->get_marks();
    }

    /// Return the mark value of dart a given mark number.
    bool get_dart_mark(const T& storage,
                        Dart_const_descriptor ADart, size_type amark) const
    {
      CGAL_assertion( ADart!=nullptr );
      return ADart->get_mark(amark);
    }

    /// Set the mark of a given mark number to a given value.
    void set_dart_mark(const T& storage,
                        Dart_const_descriptor ADart, size_type amark, bool avalue) const
    {
      CGAL_assertion( ADart!=nullptr );
      ADart->set_mark(amark, avalue);
    }

    /// Flip the mark of a given mark number to a given value.
    void flip_dart_mark(const T& storage,
                        Dart_const_descriptor ADart, size_type amark) const
    {
      CGAL_assertion( ADart!=nullptr );
      ADart->flip_mark(amark);
    }

    bool get_mask(size_type amark) const
    { return mmask_marks[amark]; }

    void flip_mask(size_type amark) const
    { mmask_marks.flip(amark); }

    /** Set simultaneously all the marks of a given dart.
     * @param adart the dart.
     * @param amarks the marks to set.
     */
    void set_marks(const T& storage,
                        Dart_const_descriptor adart,
                   const std::bitset<NB_MARKS> & amarks) const
    { set_dart_marks(storage, adart, amarks ^ mmask_marks); }

    /** Get simultaneously all the marks of a given dart.
     * @param adart the dart.
     * @return allt the marks of adart.
     */
    std::bitset<NB_MARKS> get_marks(const T& storage,
                        Dart_const_descriptor adart) const
    { return get_dart_marks(storage, adart) ^ mmask_marks; }

    /** Get the mask associated to a given mark.
     * @param amark the mark.
     * @return the mask associated to mark amark.
     */
    bool get_mask_mark(size_type amark) const
    {
      CGAL_assertion(amark>=0 && amark<NB_MARKS);
      return mmask_marks[amark];
    }

  protected:
    /// Mask marks to know the value of unmark dart, for each index i.
    mutable std::bitset<NB_MARKS> mmask_marks;
  };

  template<typename T> // TODO what is T ??
  class Mark_management_bitset_on_dart_with_index
  {
  public:
    using size_type=typename T::size_type;
    using Dart_descriptor=typename T::Dart_descriptor;
    using Dart_const_descriptor=typename T::Dart_const_descriptor;

    /// Number of marks
    static const size_type NB_MARKS = 32;

    void reset()
    { mmask_marks.reset(); }

    void swap(Mark_management_bitset_on_dart_with_index<T>& other)
    { std::swap(mmask_marks,other.mmask_marks); }

    /// init_dart_mark is called when a new dart is created.
    void init_dart_mark(const T& storage,
                        Dart_descriptor ADart)
    { set_dart_marks(storage, ADart, mmask_marks); }

    /// on_erase_dart is called when a dart is deleted
    void on_erase_dart(const T& storage,
                        Dart_descriptor ADart)
    {}

    /// Set simultaneously all the marks of this dart to a given value.
    void set_dart_marks(const T& storage,
                        Dart_const_descriptor ADart,
                        const std::bitset<NB_MARKS>& amarks) const
    { storage.mdarts[ADart].set_marks(amarks); }

    /// Return all the marks of a dart.
    std::bitset<NB_MARKS> get_dart_marks(const T& storage,
                                         Dart_const_descriptor ADart) const
    { return storage.mdarts[ADart].get_marks(); }

    /// Return the mark value of dart a given mark number.
    bool get_dart_mark(const T& storage,
                       Dart_const_descriptor ADart, size_type amark) const
    { return storage.mdarts[ADart].get_mark(amark); }

    /// Set the mark of a given mark number to a given value.
    void set_dart_mark(const T& storage,
                       Dart_const_descriptor ADart, size_type amark, bool avalue) const
    { storage.mdarts[ADart].set_mark(amark, avalue); }

    /// Flip the mark of a given mark number to a given value.
    void flip_dart_mark(const T& storage,
                        Dart_const_descriptor ADart, size_type amark) const
    { storage.mdarts[ADart].flip_mark(amark); }

    bool get_mask(size_type amark) const
    { return mmask_marks[amark]; }

    void flip_mask(size_type amark) const
    { mmask_marks.flip(amark); }

    /** Set simultaneously all the marks of a given dart.
     * @param adart the dart.
     * @param amarks the marks to set.
     */
    void set_marks(const T& storage,
                   Dart_const_descriptor adart,
                   const std::bitset<NB_MARKS> & amarks) const
    { set_dart_marks(storage, adart, amarks ^ mmask_marks); }

    /** Get simultaneously all the marks of a given dart.
     * @param adart the dart.
     * @return allt the marks of adart.
     */
    std::bitset<NB_MARKS> get_marks(const T& storage,
                                    Dart_const_descriptor adart) const
    { return get_dart_marks(storage, adart) ^ mmask_marks; }

    /** Get the mask associated to a given mark.
     * @param amark the mark.
     * @return the mask associated to mark amark.
     */
    bool get_mask_mark(size_type amark) const
    {
      CGAL_assertion(amark>=0 && amark<NB_MARKS);
      return mmask_marks[amark];
    }

  protected:
    /// Mask marks to know the value of unmark dart, for each index i.
    mutable std::bitset<NB_MARKS> mmask_marks;
  };

  template<typename T> // TODO what is T ??
  class Mark_management_bool_vector_in_map
  {
  public:

  protected:
  };

} // namespace CGAL

#endif // CGAL_MARK_MANAGEMENT_H //
// EOF //
