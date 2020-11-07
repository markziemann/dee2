module SearchPage.Types exposing (..)

import Array exposing (Array)
import Dict
import Http
import SharedTypes exposing (PaginationOffset, WebData)


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
    { searchResults: WebData SearchResults
    , searchMode : SearchMode
    , searchString : String
    , paginationOffset : SharedTypes.PaginationOffset
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
    , searchMode : SearchMode
    , searchSuggestions : WebData SearchSuggestions
    , activeSuggestion : Maybe Int
    , suggestionsVisible : Bool
    }


type Msg
    = SearchUpdate String
    | Search SearchMode String SharedTypes.PaginationOffset
    | GetSearchSuggestions String
    | GotSearchSuggestions (WebData SearchSuggestions)
    | ArrowUp
    | ArrowDown
    | EnterKey
    | StrictSelected String
    | FuzzySelected String
    | SuggestionSelected Int
    | ClickOutOfSuggestions
    | GotHttpSearchResponse SharedTypes.PaginationOffset (WebData SearchResults)
