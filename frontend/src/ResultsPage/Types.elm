module ResultsPage.Types exposing (..)

import Array
import Browser.Navigation as Nav
import Dict
import SearchPage.Types exposing (SearchParameters, SearchResults)
import Set
import SharedTypes exposing (WebData)
import Table


type alias SelectedResult =
    ( Int, ( String, String ) )


type alias SelectedResults =
    Dict.Dict Int ( String, String )


type alias ResultsPendingRemoval =
    Set.Set Int


type alias ResultRows =
    Array.Array SearchPage.Types.SearchResult

type MaybeExpired a
    = Current a
    | Expired a

type alias Model =
    { navKey : Nav.Key
    , searchResults : MaybeExpired (WebData SearchResults)
    , searchParameters: Maybe SearchParameters
    , resultsTableQuery : String
    , resultsTableState : Table.State
    , selectedResultsTableQuery : String
    , selectedResultsTableState : Table.State
    , downloading : Bool
    , selectedResults : SelectedResults
    , resultsPendingRemoval : ResultsPendingRemoval

    }


type Msg
    = ResultClicked SearchPage.Types.SearchResult
    | RemoveStagedSelections
    | SelectedResultClicked Int
    | SetResultsTableQuery String
    | SetResultsTableState Table.State
    | SetSelectedResultsTableQuery String
    | SetSelectedResultsTableState Table.State
    | DownloadRequested
    | DownloadButtonReset
    | PageRequest SearchParameters
    | ShowToggleTip Int

