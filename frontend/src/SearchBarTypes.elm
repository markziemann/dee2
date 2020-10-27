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
    , searchResults: SearchResults
    }

type Key
    = ArrowUp
    | ArrowDown

type Msg = SearchUpdate String
    | Search
    | GetSearchSuggestions String
    | GotSearchSuggestions (Result Http.Error SearchSuggestions)
    | KeyPressed Key
    | SuggestionSelected Int
    | GotHttpSearchResponse (Result Http.Error SearchResults)
    | ResultClicked (SearchResults -> SearchResults)