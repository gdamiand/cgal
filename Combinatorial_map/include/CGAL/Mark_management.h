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

  template<typename Refs>
  class Mark_management_fix
  {
  public:
    using size_type=typename Refs::size_type;

    /// Number of marks
    static const size_type NB_MARKS = 32;

    /** Count the number of used marks.
     * @return the number of used marks.
     */
    size_type number_of_used_marks() const
    { return mnb_used_marks; }

    /** Tests if a given mark is reserved.
     *  @return true iff the mark is reserved (i.e. in used).
     */
    bool is_reserved(size_type amark) const
    {
      CGAL_assertion(amark<NB_MARKS);
      return (mnb_times_reserved_marks[amark]!=0);
    }

    /**  Count the number of marked darts for a given mark.
     * @param amark the mark index.
     * @return the number of marked darts for amark.
     */
    size_type number_of_marked_darts(size_type amark) const
    {
      CGAL_assertion( is_reserved(amark) );
      return mnb_marked_darts[amark];
    }

    /**  Count the number of unmarked darts for a given mark.
     * @param amark the mark index.
     * @return the number of unmarked darts for amark.
     */
    size_type number_of_unmarked_darts(size_type amark) const
    {
      CGAL_assertion( is_reserved(amark) );
      return number_of_darts() - number_of_marked_darts(amark);
    }

    /** Tests if all the darts are unmarked for a given mark.
     * @param amark the mark index.
     * @return true iff all the darts are unmarked for amark.
     */
    bool is_whole_map_unmarked(size_type amark) const
    { return number_of_marked_darts(amark) == 0; }

    /** Tests if all the darts are marked for a given mark.
     * @param amark the mark index.
     * @return true iff all the darts are marked for amark.
     */
    bool is_whole_map_marked(size_type amark) const
    {  return number_of_marked_darts(amark) == number_of_darts(); }

    /** Reserve a new mark.
     * Get a new free mark and return its index.
     * All the darts are unmarked for this mark.
     * @return the index of the new mark.
     * @pre mnb_used_marks < NB_MARKS
     */
    size_type get_new_mark() const
    {
      if (mnb_used_marks == NB_MARKS)
      {
        std::cerr << "Not enough Boolean marks: "
                     "increase NB_MARKS in item class." << std::endl;
        std::cerr << "  (exception launched)" << std::endl;
        throw Exception_no_more_available_mark();
      }

      size_type m = mfree_marks_stack[mnb_used_marks];
      mused_marks_stack[mnb_used_marks] = m;

      mindex_marks[m] = mnb_used_marks;
      mnb_times_reserved_marks[m]=1;

      ++mnb_used_marks;
      mmark_management.on_get_new_mark(*this, m);
      CGAL_assertion(is_whole_map_unmarked(m));

      return m;
    }

    /** Increase the number of times a mark is reserved.
     *  @param amark the mark to share.
     */
    void share_a_mark(size_type amark) const
    {
      CGAL_assertion( is_reserved(amark) );
      ++mnb_times_reserved_marks[amark];
    }

    /** @return the number of times a mark is reserved.
     *  @param amark the mark to share.
     */
    size_type get_number_of_times_mark_reserved(size_type amark) const
    {
      CGAL_assertion( amark<NB_MARKS );
      return mnb_times_reserved_marks[amark];
    }

    /** Negate the mark of all the darts for a given mark.
     * After this call, all the marked darts become unmarked and all the
     * unmarked darts become marked (in constant time operation).
     * @param amark the mark index
     */
    void negate_mark(size_type amark) const
    {
      CGAL_assertion( is_reserved(amark) );

      mnb_marked_darts[amark] = number_of_darts() - mnb_marked_darts[amark];
      mmark_management.flip_mask(amark);
    }

    /** Tests if a given dart is marked for a given mark.
     * @param adart the dart to test.
     * @param amark the given mark.
     * @return true iff adart is marked for the mark amark.
     */
    bool is_marked(Dart_const_descriptor adart, size_type amark) const
    {
      CGAL_assertion( is_reserved(amark) );

      return mmark_management.get_dart_mark(*this, adart, amark)!=
             mmark_management.get_mask_mark(amark);
    }

    /** Set the mark of a given dart to a state (on or off).
     * @param adart the dart.
     * @param amark the given mark.
     * @param astate the state of the mark (on or off).
     */
    void set_mark_to(Dart_const_descriptor adart, size_type amark,
                     bool astate) const
    {
      CGAL_assertion( adart != null_dart_descriptor );
      CGAL_assertion( is_reserved(amark) );

      if (is_marked(adart, amark) != astate)
      {
        if (astate) ++mnb_marked_darts[amark];
        else --mnb_marked_darts[amark];

        flip_dart_mark(adart, amark);
      }
    }

    /** Mark the given dart.
     * @param adart the dart.
     * @param amark the given mark.
     */
    void mark(Dart_const_descriptor adart, size_type amark) const
    {
      CGAL_assertion( adart != null_dart_descriptor );
      CGAL_assertion( is_reserved(amark) );

      if (is_marked(adart, amark)) return;

      ++mnb_marked_darts[amark];
      mmark_management.flip_dart_mark(*this, adart, amark);
    }

    /** Unmark the given dart.
     * @param adart the dart.
     * @param amark the given mark.
     */
    void unmark(Dart_const_descriptor adart, size_type amark) const
    {
      CGAL_assertion( adart != null_dart_descriptor );
      CGAL_assertion( is_reserved(amark) );

      if (!is_marked(adart, amark)) return;

      --mnb_marked_darts[amark];
      mmark_management.flip_dart_mark(*this, adart, amark);
    }

    /** Mark null_dart (used as a sentinel in iterators).
     * As null dart does not belong to the set of darts, it is not counted
     * as number of marked darts.
     * @param amark the given mark.
     */
    void mark_null_dart(size_type amark) const
    {
      CGAL_assertion( is_reserved(amark) );
      mmark_management.set_dart_mark(*this, null_dart_descriptor, amark,
                                     !mmark_management.get_mask_mark(amark));
    }

    /** Unmark null_dart.
     * @param amark the given mark.
     */
    void unmark_null_dart(size_type amark) const
    {
      CGAL_assertion( is_reserved(amark) );
      mmark_management.set_dart_mark(*this, null_dart_descriptor, amark,
                                     mmark_management.get_mask_mark(amark));
    }

    /** Unmark all the darts of the map for a given mark.
     * If all the darts are marked or unmarked, this operation takes \cgalBigO{1}
     * operations, otherwise it traverses all the darts of the map.
     * @param amark the given mark.
     */
    void unmark_all(size_type amark) const
    {
      CGAL_assertion( is_reserved(amark) );

      if ( is_whole_map_marked(amark) )
      {
        negate_mark(amark);
      }
      else if ( !is_whole_map_unmarked(amark) )
      {
        for ( typename Dart_range::const_iterator it(darts().begin()),
             itend(darts().end()); it!=itend; ++it)
          unmark(it, amark);
      }
      CGAL_assertion(is_whole_map_unmarked(amark));
      unmark_null_dart(amark);
    }

    /** Free a given mark, previously calling unmark_all_darts.
     * @param amark the given mark.
     */
    void free_mark(size_type amark) const
    {
      CGAL_assertion( is_reserved(amark) );

      if ( mnb_times_reserved_marks[amark]>1 )
      {
        --mnb_times_reserved_marks[amark];
        return;
      }

      unmark_all(amark);

      // 1) We remove amark from the array mused_marks_stack by
      //    replacing it with the last mark in this array.
      mused_marks_stack[mindex_marks[amark]] =
          mused_marks_stack[--mnb_used_marks];
      mindex_marks[mused_marks_stack[mnb_used_marks]] =
          mindex_marks[amark];

      // 2) We add amark in the array mfree_marks_stack and update its index.
      mfree_marks_stack[ mnb_used_marks ] = amark;
      mindex_marks[amark] = mnb_used_marks;

      mnb_times_reserved_marks[amark]=0;

      mmark_management.on_free_mark(*this, amark);
    }


  protected:
    /// Number of times each mark is reserved. 0 if the mark is free.
    mutable size_type mnb_times_reserved_marks[NB_MARKS];

    /// Number of used marks.
    mutable size_type mnb_used_marks;

    /// Index of each mark, in mfree_marks_stack or in mfree_marks_stack.
    mutable size_type mindex_marks[NB_MARKS];

    /// "Stack" of free marks.
    mutable size_type mfree_marks_stack[NB_MARKS];

    /// "Stack" of used marks.
    mutable size_type mused_marks_stack[NB_MARKS];

    /// Number of marked darts for each used marks.
    mutable size_type mnb_marked_darts[NB_MARKS];
  };

  template<typename Refs>
  class Mark_management_bitset_on_dart
  {
  public:
    using size_type=typename Refs::size_type;
    using Dart_descriptor=typename Refs::Dart_descriptor;
    using Dart_const_descriptor=typename Refs::Dart_const_descriptor;


    void reset()
    { mmask_marks.reset(); }

    void swap(Mark_management_bitset_on_dart<Refs>& other)
    { std::swap(mmask_marks, other.mmask_marks); }

    /// on_new_dart is called when a new dart is created.
    void on_new_dart(const Refs& storage, Dart_descriptor ADart) const
    { set_dart_marks(storage, ADart, mmask_marks); }

    /// on_delete_dart is called when a dart is deleted
    void on_delete_dart(const Refs& storage, Dart_descriptor ADart) const
    {}

    /// on_get_new_mark is called when a new mark is reserved.
    void on_get_new_mark(const Refs& storage, size_type amark) const
    {}

    /// on_free_mark is called when a new mark is released.
    void on_free_mark(const Refs& storage, size_type amark) const
    {}

    /// Return the mark value of dart a given mark number.
    bool get_dart_mark(const Refs& storage,
                       Dart_const_descriptor ADart, size_type amark) const
    {
      CGAL_assertion( ADart!=nullptr );
      return ADart->get_mark(amark);
    }

    /// Set the mark of a given mark number to a given value.
    void set_dart_mark(const Refs& storage,
                       Dart_const_descriptor ADart,
                       size_type amark, bool avalue) const
    {
      CGAL_assertion( ADart!=nullptr );
      ADart->set_mark(amark, avalue);
    }

    /// Flip the mark of a given mark number to a given value.
    void flip_dart_mark(const Refs& storage,
                        Dart_const_descriptor ADart, size_type amark) const
    {
      CGAL_assertion( ADart!=nullptr );
      ADart->flip_mark(amark);
    }

    void flip_mask(size_type amark) const
    { mmask_marks.flip(amark); }

    /// Set simultaneously all the marks of this dart to a given value.
    void set_dart_marks(const Refs& storage,
                        Dart_const_descriptor ADart,
                        const std::bitset<NB_MARKS>& amarks) const
    {
      CGAL_assertion( ADart!=nullptr );
      ADart->set_marks(amarks);
    }

    /// Return all the marks of a dart.
    std::bitset<NB_MARKS> get_dart_marks(const Refs& storage,
                                         Dart_const_descriptor ADart) const
    {
      CGAL_assertion( ADart!=nullptr );
      return ADart->get_marks();
    }

    /** Set simultaneously all the marks of a given dart.
     * @param adart the dart.
     * @param amarks the marks to set.
     */
    void set_marks(const Refs& storage,
                   Dart_const_descriptor adart,
                   const std::bitset<NB_MARKS> & amarks) const
    { set_dart_marks(storage, adart, amarks ^ mmask_marks); }

    /** Get simultaneously all the marks of a given dart.
     * @param adart the dart.
     * @return allt the marks of adart.
     */
    std::bitset<NB_MARKS> get_marks(const Refs& storage,
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

  template<typename Refs>
  class Mark_management_bitset_on_dart_with_index
  {
  public:
    using size_type=typename Refs::size_type;
    using Dart_descriptor=typename Refs::Dart_descriptor;
    using Dart_const_descriptor=typename Refs::Dart_const_descriptor;

    /// Number of marks
    static const size_type NB_MARKS = 32;

    void reset()
    { mmask_marks.reset(); }

    void swap(Mark_management_bitset_on_dart_with_index<Refs>& other)
    { std::swap(mmask_marks, other.mmask_marks); }

    /// on_new_dart is called when a new dart is created.
    void on_new_dart(const Refs& storage, Dart_descriptor ADart) const
    { set_dart_marks(storage, ADart, mmask_marks); }

    /// on_delete_dart is called when a dart is deleted
    void on_delete_dart(const Refs& storage, Dart_descriptor ADart) const
    {}

    /// on_get_new_mark is called when a new mark is reserved.
    void on_get_new_mark(const Refs& storage, size_type amark) const
    {}

    /// on_free_mark is called when a new mark is released.
    void on_free_mark(const Refs& storage, size_type amark) const
    {}

    /// Return the mark value of dart a given mark number.
    bool get_dart_mark(const Refs& storage,
                       Dart_const_descriptor ADart, size_type amark) const
    { return storage.mdarts[ADart].get_mark(amark); }

    /// Set the mark of a given mark number to a given value.
    void set_dart_mark(const Refs& storage,
                       Dart_const_descriptor ADart,
                       size_type amark, bool avalue) const
    { storage.mdarts[ADart].set_mark(amark, avalue); }

    /// Flip the mark of a given mark number to a given value.
    void flip_dart_mark(const Refs& storage,
                        Dart_const_descriptor ADart, size_type amark) const
    { storage.mdarts[ADart].flip_mark(amark); }

    void flip_mask(size_type amark) const
    { mmask_marks.flip(amark); }

    /// Set simultaneously all the marks of this dart to a given value.
    void set_dart_marks(const Refs& storage,
                        Dart_const_descriptor ADart,
                        const std::bitset<NB_MARKS>& amarks) const
    { storage.mdarts[ADart].set_marks(amarks); }

    /// Return all the marks of a dart.
    std::bitset<NB_MARKS> get_dart_marks(const Refs& storage,
                                         Dart_const_descriptor ADart) const
    { return storage.mdarts[ADart].get_marks(); }

    /** Set simultaneously all the marks of a given dart.
     * @param adart the dart.
     * @param amarks the marks to set.
     */
    void set_marks(const Refs& storage,
                   Dart_const_descriptor adart,
                   const std::bitset<NB_MARKS> & amarks) const
    { set_dart_marks(storage, adart, amarks ^ mmask_marks); }

    /** Get simultaneously all the marks of a given dart.
     * @param adart the dart.
     * @return allt the marks of adart.
     */
    std::bitset<NB_MARKS> get_marks(const Refs& storage,
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

  template<typename Refs>
  class Mark_management_bool_vector_in_map
  {
  public:
    using size_type=typename Refs::size_type;
    using Dart_descriptor=typename Refs::Dart_descriptor;
    using Dart_const_descriptor=typename Refs::Dart_const_descriptor;

    /// Number of marks
    static const size_type NB_MARKS = 32;

    void reset()
    {
      mmask_marks.reset();
      for(auto& v: marray_of_marks)
      { v.clear(); }
    }

    void swap(Mark_management_bool_vector_in_map<Refs>& other)
    {
      std::swap(mmask_marks, other.mmask_marks);
      std::swap(marray_of_marks, other.marray_of_marks);
    }

    /// on_new_dart is called when a new dart is created.
    void on_new_dart(const Refs& storage, Dart_descriptor ADart)
    {}

    /// on_delete_dart is called when a dart is deleted
    void on_delete_dart(const Refs& storage, Dart_descriptor ADart)
    {}

    /// on_get_new_mark is called when a new mark is reserved.
    void on_get_new_mark(const Refs& storage, size_type amark) const
    {}

    /// on_free_mark is called when a new mark is released.
    void on_free_mark(const Refs& storage, size_type amark) const
    {}

    /// Return the mark value of dart a given mark number.
    bool get_dart_mark(const Refs& storage,
                       Dart_const_descriptor ADart, size_type amark) const
    {}

    /// Set the mark of a given mark number to a given value.
    void set_dart_mark(const Refs& storage,
                       Dart_const_descriptor ADart,
                       size_type amark, bool avalue) const
    {}

    /// Flip the mark of a given mark number to a given value.
    void flip_dart_mark(const Refs& storage,
                        Dart_const_descriptor ADart, size_type amark) const
    {}

    void flip_mask(size_type amark) const
    { mmask_marks.flip(amark); }

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

    /// array of vector of marks
    std::array<std::vector<bool>, NB_MARKS> marray_of_marks;
  };

} // namespace CGAL

#endif // CGAL_MARK_MANAGEMENT_H //
// EOF //
