// $Id$
//==============================================================================
//!
//! \file XMLInputBase.h
//!
//! \date Jul 16 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Base class for XML input parsing functionality.
//!
//==============================================================================

#ifndef _XML_INPUT_BASE_H
#define _XML_INPUT_BASE_H

namespace tinyxml2 {
  class XMLDocument;
  class XMLElement;
}


/*!
  \brief Base class for XML-based input file parsing.
  \details This class is inherited by SIMadmin for input parsing handling,
  and can also be used by applications for pre-parsing of the input file.
*/

class XMLInputBase
{
public:
  //! \brief Reads an XML input file.
  //! \param[in] fName Name of the file to read
  //! \param[in] verbose If \e true, print the tags being parsed
  bool readXML(const char* fName, bool verbose = true);

  //! \brief Loads data from an XML-formatted text string.
  //! \details This method is a convenience offered for unit testing only.
  bool loadXML(const char* xml);

protected:
  //! \brief Loads an XML-file into a tinyxml2::XMLDocument object.
  //! \param[out] doc The top-level node of the XML-file
  //! \param[in] fName Name of the file to read
  //! \param[in] verbose If \e true, print the file(s) being read
  //! \return Pointer to the root element of the XML-file
  const tinyxml2::XMLElement* loadFile(tinyxml2::XMLDocument& doc,
                                       const char* fName, bool verbose = false);

  //! \brief Parses a data section from an XML element.
  virtual bool parse(const tinyxml2::XMLElement* elem) = 0;

  //! \brief Returns a list of prioritized XML-tags.
  //! \details Certain tags need to be parsed before others.
  //! This method will take care of this. It is called by the readXML() method
  //! in order to read some top level tags in a certain order.
  virtual const char** getPrioritizedTags() const { return nullptr; }
};

#endif
