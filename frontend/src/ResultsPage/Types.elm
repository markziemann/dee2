module ResultsPage.Types exposing (..)

import Array
import SearchPage.Types
import Table


type alias Model =
    { searchHits : Maybe Int
    , searchResultRows : Maybe (Array.Array SearchPage.Types.SearchResult)
    , resultsTableState : Table.State
    , resultsTableQuery : String
    , downloading : Bool
    }


type Msg
    = ResultClicked SearchPage.Types.SearchResult
    | SetResultsTableQuery String
    | SetResultsTableState Table.State
    | DownloadRequested
    | DownloadButtonReset
