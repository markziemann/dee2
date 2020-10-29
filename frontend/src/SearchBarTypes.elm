module SearchBarTypes exposing (..)

import Array exposing (Array)
import Http


type alias SearchData =
    List ( String, String )


type alias SearchResult =
    { id : Int
    , data : SearchData
    , selected : Bool
    }


type alias SearchResults =
    Array SearchResult


type alias SearchSuggestions =
    Array String


type alias ActiveSuggestion =
    Int


type alias Model =
    { searchString : String
    , searchSuggestions : SearchSuggestions
    , activeSuggestion : Maybe Int
    , suggestionsVisible : Bool
    , waitingForResponse : Bool
    }


type Msg
    = SearchUpdate String
    | Search
    | GetSearchSuggestions String
    | GotSearchSuggestions (Result Http.Error SearchSuggestions)
    | ArrowUp
    | ArrowDown
    | SuggestionSelected Int
    | ClickOutOfSuggestions
    | GotHttpSearchResponse (Result Http.Error SearchResults)
