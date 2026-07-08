// $Id$
//==============================================================================
//!
//! \file XMLInputBase.C
//!
//! \date Jul 16 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Base class for XML input parsing functionality.
//!
//==============================================================================

#include "XMLInputBase.h"
#include "IFEM.h"
#include "tinyxml2.h"

#include <algorithm>
#include <cstring>
#include <iomanip>
#include <vector>


namespace
{
  //! \brief Helper class to handle included files.
  class IncludeInjector : public tinyxml2::XMLVisitor
  {
  public:
    //! \brief The constructor creates the new document.
    IncludeInjector(int i) : new_doc(true,tinyxml2::COLLAPSE_WHITESPACE)
    {
      level = i;
    }

    //! \brief Call-back invoked when visiting a new xml-tag.
    //! \details If \a elem is an &lt;include&gt; tag,
    //! it is replaced by the content of the included XML-file,
    //! otherwise its content is copied to \ref new_doc.
    bool VisitEnter(const tinyxml2::XMLElement& elem,
                    const tinyxml2::XMLAttribute* attribute) override
    {
      const char* str = elem.GetText();
      if (!strcmp(elem.Value(),"include") && str)
      {
        tinyxml2::XMLDocument doc(true,tinyxml2::COLLAPSE_WHITESPACE);
        if (doc.LoadFile(str) != tinyxml2::XML_SUCCESS)
        {
          std::cerr <<"** XMLInputBase: Failed to load included file \""
                    << str <<"\".\n"<< doc.ErrorStr() << std::endl;
          return false;
        }

        if (level > 0)
          IFEM::cout << std::setw(2*level+7) <<"Loading"
                     <<" included file "<< str << std::endl;
        for (const tinyxml2::XMLElement* child = doc.FirstChildElement();
             child; child = child->NextSiblingElement())
          currElem->InsertEndChild(child->DeepClone(&new_doc));

        include_found = include_processed = true;
      }
      else
      {
        tinyxml2::XMLElement* e = new_doc.NewElement(elem.Name());
        if (str)
          e->SetText(str);
        if (!currElem)
          new_doc.InsertEndChild(e);
        else
          currElem->InsertEndChild(e);
        while (attribute)
        {
          e->SetAttribute(attribute->Name(),attribute->Value());
          attribute = attribute->Next();
        }
        currElem = e;
      }
      return true;
    }

    //! \brief Call-back invoked when leaving an xml-tag.
    //! \details This method sets the element to insert into to the parent.
    bool VisitExit(const tinyxml2::XMLElement&) override
    {
      if (include_processed)
        include_processed = false;
      else if (currElem)
        if (tinyxml2::XMLNode* parent = currElem->Parent(); parent)
          currElem = parent->ToElement();
      return true;
    }

    //! \brief Inserts content from included file (if any) into \a doc.
    bool replaceIncluded(tinyxml2::XMLDocument& doc)
    {
      if (!include_found)
        return false;

      new_doc.DeepCopy(&doc);
      return true;
    }

  private:
    bool include_found     = false; //!< If \e true, an include tag is found
    bool include_processed = false; //!< If \e true, an include tag is processed
    tinyxml2::XMLElement* currElem = nullptr; //!< Current root element
    tinyxml2::XMLDocument new_doc; //!< XML content with include file replaced
    int level; //!< For indented print of file names
  };
}


const tinyxml2::XMLElement* XMLInputBase::loadFile (tinyxml2::XMLDocument& doc,
                                                    const char* fName,
                                                    bool verbose)
{
  if (doc.LoadFile(fName) != tinyxml2::XML_SUCCESS)
  {
    std::cerr <<" *** XMLInputBase: Failed to load XML-file \""<< fName
              <<"\".\n" << doc.ErrorStr() << std::endl;
    return nullptr;
  }

  const tinyxml2::XMLElement* tag = doc.RootElement();
  if (!tag || !tag->Value())
    return nullptr; // Empty file?

  if (strcmp(tag->Value(),"simulation"))
  {
    std::cerr <<" *** XMLInputBase: Malformatted XML-file \""<< fName
              <<"\".\n     The root tag must be <simulation> - not <"
              << tag->Value() <<">"<< std::endl;
    return nullptr;
  }

  if (verbose)
    IFEM::cout <<"\nParsing input file "<< fName << std::endl;

  for (size_t i = 0; i < 10 && tag; i++) // Maximum 10 levels of include files
    if (IncludeInjector v(verbose ? i+1 : 0); !doc.RootElement()->Accept(&v))
      tag = nullptr; // abort on include failure
    else if (!v.replaceIncluded(doc))
      break; // no include tags

#ifdef SP_DEBUG
  if (verbose) {
    std::cout <<"\nHere is the input-file content:"<< std::endl;
    doc.Print();
  }
#endif

  return tag ? doc.RootElement() : nullptr;
}


bool XMLInputBase::readXML (const char* fName, bool verbose)
{
  tinyxml2::XMLDocument doc(true, tinyxml2::COLLAPSE_WHITESPACE);
  const tinyxml2::XMLElement* tag = this->loadFile(doc,fName,verbose);
  if (!tag) return false;

  // Lambda function for parsing an XML-tag with logging and failure message.
  auto parseTag = [this,verbose](const tinyxml2::XMLElement* tag)
  {
    if (verbose)
      IFEM::cout <<"\nParsing <"<< tag->Value() <<">"<< std::endl;

    if (this->parse(tag))
      return true;

    std::cerr <<" *** XMLInputBase: Failure occurred while parsing <"
              << tag->Value() <<">."<< std::endl;
    return false;
  };

  std::vector<const tinyxml2::XMLElement*> parsed;
  if (const char** q = this->getPrioritizedTags(); q)
    while (*q)
      if (const tinyxml2::XMLElement* elm = tag->FirstChildElement(*(q++)); elm)
      {
        if (!parseTag(elm))
          return false;
        parsed.push_back(elm);
      }

  for (tag = tag->FirstChildElement(); tag; tag = tag->NextSiblingElement())
    if (std::find(parsed.begin(),parsed.end(),tag) == parsed.end())
      if (!parseTag(tag))
        return false;

  if (verbose)
    IFEM::cout <<"\nParsing input file succeeded."<< std::endl;

  return true;
}


bool XMLInputBase::loadXML (const char* xml)
{
  tinyxml2::XMLDocument doc;
  doc.Parse(xml);
  const tinyxml2::XMLElement* tag = doc.RootElement();
  return tag ? this->parse(tag) : false;
}
