module SearchPage.Types exposing (..)

import Array exposing (Array)
import Browser.Navigation as Nav
import Dict
import SharedTypes exposing (PaginationOffset, WebData)

type SearchParameters
    = SearchParameters SearchMode String SharedTypes.PaginationOffset

type alias SearchData =
    Dict.Dict String String


type alias SearchResult =
    { id : Int
    , data : SearchData
    }


type alias SearchResults =
    { hits : Int
    , rows : Array SearchResult
    }





type alias SearchSuggestions =
    Array String


type alias ActiveSuggestion =
    Int


type SearchMode
    = Strict
    | Fuzzy


type alias Model =
    { navKey : Nav.Key
    , searchParameters : SearchParameters
    , defaultPaginationOffset: PaginationOffset
    , searchSuggestions : WebData SearchSuggestions
    , activeSuggestion : Maybe Int
    , suggestionsVisible : Bool
    }


type Msg
    = SearchUpdate String
    | Search SearchParameters
    | GetSearchSuggestions SearchParameters
    | GotSearchSuggestions (WebData SearchSuggestions)
    | ArrowUp
    | ArrowDown
    | EnterKey
    | StrictSelected String
    | FuzzySelected String
    | SuggestionSelected Int
    | ClickOutOfSuggestions
