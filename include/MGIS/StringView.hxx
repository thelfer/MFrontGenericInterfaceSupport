/*!
 * \file   StringView.hxx
 * \brief
 * \author Thomas Helfer
 * \date   02/08/2018
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 *
 * Based on the work of Matthew Rodusek (matthew.rodusek@gmail.com)
 * distributed under the MIT License (MIT). See
 * https://github.com/bitwizeshift/string_view-standalone for details.
 *
 * Licensed under the MIT License <http://opensource.org/licenses/MIT>.
 * Copyright (c) 2016 Matthew Rodusek <http://rodusek.me>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef LIB_MGIS_STRINGVIEW_HXX
#define LIB_MGIS_STRINGVIEW_HXX

#include <algorithm>
#include <string>
#include <ostream>

namespace mgis{

  ////////////////////////////////////////////////////////////////////////////
  /// \brief A wrapper around non-owned strings.
  ///
  /// This is an implementation of the C++17 string_view proposal
  ////////////////////////////////////////////////////////////////////////////
  template <typename CharT, typename Traits = std::char_traits<CharT>>
  class basic_string_view final {
    //------------------------------------------------------------------------
    // Public Member Types
    //------------------------------------------------------------------------
   public:
    using char_type = CharT;
    using traits_type = Traits;
    using size_type = size_t;

    using value_type = CharT;
    using reference = value_type&;
    using const_reference = const value_type&;
    using pointer = value_type*;
    using const_pointer = const value_type*;

    using iterator = const CharT*;
    using const_iterator = const CharT*;

    //------------------------------------------------------------------------
    // Public Members
    //------------------------------------------------------------------------
   public:
    static constexpr size_type npos = size_type(-1);

    //------------------------------------------------------------------------
    // Constructors
    //------------------------------------------------------------------------
   public:
    /// \brief Default constructs a basic_string_view without any content
    constexpr basic_string_view() noexcept;

    /// \brief Constructs a basic_string_view by copying another one
    ///
    /// \param other the string view being copied
    constexpr basic_string_view(const basic_string_view& other) noexcept =
        default;

    /// \brief Constructs a basic_string_view by moving anothe rone
    ///
    /// \param other the string view being moved
    constexpr basic_string_view(basic_string_view&& other) noexcept = default;

    /// \brief Constructs a basic_string_view from a std::basic_string
    ///
    /// \param str the string to view
    template <typename Allocator>
    constexpr basic_string_view(
        const std::basic_string<CharT, Traits, Allocator>& str) noexcept;

    /// \brief Constructs a basic_string_view from an ansi-string
    ///
    /// \param str the string to view
    constexpr basic_string_view(const char_type* str) noexcept;

    /// \brief Constructs a basic_string_view from an ansi string of a given
    /// size
    ///
    /// \param str the string to view
    /// \param count the size of the string
    constexpr basic_string_view(const char_type* str, size_type count) noexcept;

    //------------------------------------------------------------------------
    // Assignment
    //------------------------------------------------------------------------
   public:
    /// \brief Assigns a basic_string_view from an ansi-string
    ///
    /// \param view the string to view
    /// \return reference to \c (*this)
    basic_string_view& operator=(const basic_string_view& view) = default;

    //------------------------------------------------------------------------
    // Capacity
    //------------------------------------------------------------------------
   public:
    /// \brief Returns the length of the string, in terms of bytes
    ///
    /// \return the length of the string, in terms of bytes
    constexpr size_type size() const noexcept;

    /// \copydoc basic_string_view::size
    constexpr size_type length() const noexcept;

    /// \brief The largest possible number of char-like objects that can be
    ///        referred to by a basic_string_view.
    /// \return Maximum number of characters
    constexpr size_type max_size() const noexcept;

    /// \brief Returns whether the basic_string_view is empty
    ///        (i.e. whether its length is 0).
    ///
    /// \return whether the basic_string_view is empty
    constexpr bool empty() const noexcept;

    //------------------------------------------------------------------------
    // Element Access
    //------------------------------------------------------------------------
   public:
    /// \brief Gets the ansi-string of the current basic_string_view
    ///
    /// \return the ansi-string pointer
    constexpr const char_type* c_str() const noexcept;

    /// \brief Gets the data of the current basic_string_view
    ///
    /// \note This is an alias of #c_str
    ///
    /// \return the data this basic_string_view contains
    constexpr const char_type* data() const noexcept;

    /// \brief Accesses the element at index \p pos
    ///
    /// \param pos the index to access
    /// \return const reference to the character
    constexpr const_reference operator[](size_t pos) const noexcept;

    /// \brief Accesses the element at index \p pos
    ///
    /// \param pos the index to access
    /// \return const reference to the character
    constexpr const_reference at(size_t pos) const;

    /// \brief Access the first character of the string
    ///
    /// \note Undefined behavior if basic_string_view is empty
    ///
    /// \return reference to the first character of the string
    constexpr const_reference front() const noexcept;

    /// \brief References the last character of the string
    ///
    /// \note Undefined behavior if basic_string_view is empty
    ///
    /// \return reference to the last character of the string
    constexpr const_reference back() const noexcept;

    //------------------------------------------------------------------------
    // Modifiers
    //------------------------------------------------------------------------
   public:
    /// \brief Moves the start of the view forward by n characters.
    ///
    /// The behavior is undefined if n > size().
    ///
    /// \param n number of characters to remove from the start of the view
    void remove_prefix(size_type n) noexcept;

    /// \brief Moves the end of the view back by n characters.
    ///
    /// The behavior is undefined if n > size().
    ///
    /// \param n number of characters to remove from the end of the view
    void remove_suffix(size_type n) noexcept;

    /// \brief Exchanges the view with that of v.
    ///
    /// \param v view to swap with
    void swap(basic_string_view& v) noexcept;

    //------------------------------------------------------------------------
    // Conversions
    //------------------------------------------------------------------------
   public:
    /// \brief Creates a basic_string with a copy of the content of the current
    /// view.
    ///
    /// \tparam Allocator type used to allocate internal storage
    /// \param a Allocator instance to use for allocating the new string
    ///
    /// \return A basic_string containing a copy of the characters of the
    /// current view.
    template <class Allocator = std::allocator<CharT>>
    constexpr std::basic_string<CharT, Traits, Allocator> to_string(
        const Allocator& a = Allocator()) const;

    /// \copydoc basic_string_view::to_string
    template <class Allocator>
    explicit constexpr operator std::basic_string<CharT, Traits, Allocator>()
        const;

    //------------------------------------------------------------------------
    // Operations
    //------------------------------------------------------------------------
   public:
    /// \brief Copies the substring [pos, pos + rcount) to the character string
    /// pointed
    ///        to by dest, where rcount is the smaller of count and size() -
    ///        pos.
    ///
    /// \param dest pointer to the destination character string
    /// \param count requested substring length
    /// \param pos position of the first character
    size_type copy(char_type* dest,
                   size_type count = npos,
                   size_type pos = 0) const;

    /// \brief Returns a substring of this viewed string
    ///
    /// \param pos the position of the first character in the substring
    /// \param len the length of the substring
    /// \return the created substring
    basic_string_view substr(size_t pos = 0, size_t len = npos) const;

    //------------------------------------------------------------------------

    /// \brief Compares two character sequences
    ///
    /// \param v view to compare
    /// \return negative value if this view is less than the other character
    ///         sequence, zero if the both character sequences are equal,
    ///         positive
    ///         value if this view is greater than the other character sequence.
    int compare(basic_string_view v) const noexcept;

    /// \brief Compares two character sequences
    ///
    /// \param pos   position of the first character in this view to compare
    /// \param count number of characters of this view to compare
    /// \param v     view to compare
    /// \return negative value if this view is less than the other character
    ///         sequence, zero if the both character sequences are equal,
    ///         positive
    ///         value if this view is greater than the other character sequence.
    int compare(size_type pos, size_type count, basic_string_view v) const;

    /// \brief Compares two character sequences
    ///
    /// \param pos1   position of the first character in this view to compare
    /// \param count1 number of characters of this view to compare
    /// \param v      view to compare
    /// \param pos2   position of the second character in this view to compare
    /// \param count2 number of characters of the given view to compare
    /// \return negative value if this view is less than the other character
    ///         sequence, zero if the both character sequences are equal,
    ///         positive
    ///         value if this view is greater than the other character sequence.
    int compare(size_type pos1,
                size_type count1,
                basic_string_view v,
                size_type pos2,
                size_type count2) const;

    /// \brief Compares two character sequences
    ///
    /// \param s pointer to the character string to compare to
    /// \return negative value if this view is less than the other character
    ///         sequence, zero if the both character sequences are equal,
    ///         positive
    ///         value if this view is greater than the other character sequence.
    int compare(const char_type* s) const;

    /// \brief Compares two character sequences
    ///
    /// \param pos   position of the first character in this view to compare
    /// \param count number of characters of this view to compare
    /// \param s pointer to the character string to compare to
    /// \return negative value if this view is less than the other character
    ///         sequence, zero if the both character sequences are equal,
    ///         positive
    ///         value if this view is greater than the other character sequence.
    int compare(size_type pos, size_type count, const char_type* s) const;

    /// \brief Compares two character sequences
    ///
    /// \param pos   position of the first character in this view to compare
    /// \param count1 number of characters of this view to compare
    /// \param s pointer to the character string to compare to
    /// \param count2 number of characters of the given view to compare
    /// \return negative value if this view is less than the other character
    ///         sequence, zero if the both character sequences are equal,
    ///         positive
    ///         value if this view is greater than the other character sequence.
    int compare(size_type pos,
                size_type count1,
                const char_type* s,
                size_type count2) const;

    //------------------------------------------------------------------------

    size_type find(basic_string_view v, size_type pos = 0) const;

    size_type find(char_type c, size_type pos = 0) const;

    size_type find(const char_type* s, size_type pos, size_type count) const;

    size_type find(const char_type* s, size_type pos = 0) const;

    //------------------------------------------------------------------------

    size_type rfind(basic_string_view v, size_type pos = npos) const;

    size_type rfind(char_type c, size_type pos = npos) const;

    size_type rfind(const char_type* s, size_type pos, size_type count) const;

    size_type rfind(const char_type* s, size_type pos = npos) const;

    //------------------------------------------------------------------------

    size_type find_first_of(basic_string_view v, size_type pos = 0) const;

    size_type find_first_of(char_type c, size_type pos = 0) const;

    size_type find_first_of(const char_type* s,
                            size_type pos,
                            size_type count) const;

    size_type find_first_of(const char_type* s, size_type pos = 0) const;

    //------------------------------------------------------------------------
    // Iterators
    //------------------------------------------------------------------------
   public:
    /// \brief Retrieves the begin iterator for this basic_string_view
    ///
    /// \return the begin iterator
    const_iterator begin() const noexcept;

    /// \brief Retrieves the end iterator for this basic_string_view
    ///
    /// \return the end iterator
    const_iterator end() const noexcept;

    /// \copydoc basic_string_view::begin()
    const_iterator cbegin() const noexcept;

    /// \copydoc basic_string_view::end()
    const_iterator cend() const noexcept;

    //------------------------------------------------------------------------
    // Private Member
    //------------------------------------------------------------------------
   private:
    const char_type* m_str;  ///< The internal string type
    size_type m_size;        ///< The size of this string
  };

  //--------------------------------------------------------------------------
  // Public Functions
  //--------------------------------------------------------------------------

  /// \brief Overload for ostream output of basic_string_view
  ///
  /// \param o   The output stream to print to
  /// \param str the string to print
  /// \return reference to the output stream
  template <typename CharT, typename Traits>
  std::basic_ostream<CharT, Traits>& operator<<(
      std::basic_ostream<CharT, Traits>& o,
      const basic_string_view<CharT, Traits>& str);

  template <typename CharT, typename Traits>
  void swap(basic_string_view<CharT, Traits>& lhs,
            basic_string_view<CharT, Traits>& rhs) noexcept;

  //--------------------------------------------------------------------------
  // Comparison Functions
  //--------------------------------------------------------------------------

  /// \brief Compares equality between the two basic_string_views
  ///
  /// \param lhs the left string view
  /// \param rhs the right string view
  /// \return true if the two strings are the same
  template <typename CharT, typename Traits>
  bool operator==(const basic_string_view<CharT, Traits>& lhs,
                  const basic_string_view<CharT, Traits>& rhs) noexcept;

  /// \brief Compares inequality between the two basic_string_views
  ///
  /// \param lhs the left string view
  /// \param rhs the right string view
  /// \return true if the two strings are different
  template <typename CharT, typename Traits>
  bool operator!=(const basic_string_view<CharT, Traits>& lhs,
                  const basic_string_view<CharT, Traits>& rhs) noexcept;

  /// \brief Checks if the left string is less than the right substring
  ///
  /// \param lhs the left string view
  /// \param rhs the right string view
  /// \return true if the left string has a character less than the right
  ///         string, or if the right string is shorter than the left string
  template <typename CharT, typename Traits>
  bool operator<(const basic_string_view<CharT, Traits>& lhs,
                 const basic_string_view<CharT, Traits>& rhs) noexcept;

  ///
  /// \param lhs
  /// \param rhs
  /// \return
  template <typename CharT, typename Traits>
  bool operator>(const basic_string_view<CharT, Traits>& lhs,
                 const basic_string_view<CharT, Traits>& rhs) noexcept;

  ///
  /// \param lhs
  /// \param rhs
  /// \return
  template <typename CharT, typename Traits>
  bool operator<=(const basic_string_view<CharT, Traits>& lhs,
                  const basic_string_view<CharT, Traits>& rhs) noexcept;

  ///
  /// \param lhs
  /// \param rhs
  /// \return
  template <typename CharT, typename Traits>
  bool operator>=(const basic_string_view<CharT, Traits>& lhs,
                  const basic_string_view<CharT, Traits>& rhs) noexcept;

  //--------------------------------------------------------------------------
  // Type Aliases
  //--------------------------------------------------------------------------

  using string_view = basic_string_view<char>;
  using wstring_view = basic_string_view<wchar_t>;
  using u16string_view = basic_string_view<char16_t>;
  using u32string_view = basic_string_view<char32_t>;

}  // end of namespace mgis

#include "MGIS/StringView.ixx"

#endif /* LIB_MGIS_STRINGVIEW_HXX */
