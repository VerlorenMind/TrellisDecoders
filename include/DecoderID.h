#ifndef BEAST_DECODERID_H
#define BEAST_DECODERID_H

#include <string>

enum class DecoderID
{
   ERROR,
   BEAST,
   BSD,
   KTKL,
};

static std::string idToString(DecoderID id)
{
   switch(id)
   {
       case DecoderID::BEAST:
           return "beast";
       case DecoderID::BSD:
           return "bsd";
       case DecoderID::KTKL:
           return "ktkl";
       default:
           return "";
   }
}
static DecoderID stringToId(const std::string &name)
{
   if(name == "beast")
   {
       return DecoderID::BEAST;
   }
   else if (name == "bsd")
   {
       return DecoderID::BSD;
   }
   else if (name == "ktkl")
   {
       return DecoderID::KTKL;
   }
   else
   {
       return DecoderID::ERROR;
   }
}
#endif //BEAST_DECODERID_H
