module SearchBarTypes exposing (..)

import Array exposing (Array)
import Dict
import Http


type alias SearchData =
    Dict.Dict String String


type alias SearchResult =
    { id : Int
    , data : SearchData
    , selected : Bool
    }


type alias SearchResults =
    { hits : Int
    , rows : Array SearchResult
    }

type alias OutMsg =
    { hits : Int
    , rows : Array SearchResult
    , searchString : String
    }

type alias SearchSuggestions =
    Array String


type alias ActiveSuggestion =
    Int

type SearchMode
    = Strict
    | Fuzzy

type alias Model =
    { searchString : String
    , searchMode: SearchMode
    , searchSuggestions : SearchSuggestions
    , activeSuggestion : Maybe Int
    , suggestionsVisible : Bool
    , waitingForResponse : Bool
    , test : Bool
    }


type Msg
    = SearchUpdate String
    | Search
    | GetSearchSuggestions String
    | GotSearchSuggestions (Result Http.Error SearchSuggestions)
    | ArrowUp
    | ArrowDown
    | EnterKey
    | StrictSelected String
    | FuzzySelected String
    | SuggestionSelected Int
    | ClickOutOfSuggestions
    | GotHttpSearchResponse (Result Http.Error SearchResults)
    | Test
